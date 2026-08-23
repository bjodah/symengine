# Thread scaling and the singleton reference counts

Branch `fix/singleton-false-sharing`. Written against ledger row **B-4** in
`symengine.rubi/UPSTREAMABLE.md`, which claimed that multi-threaded scaling
caps at 2.02× because "SymEngine's arithmetic singletons are process-global
and their intrusive reference counts false-share three heap cache lines".

The claim is right about *where* and wrong about *what*: the three lines are
`zero`, `one` and `minus_one`, and the traffic on each is not false sharing
between neighbours but **true sharing of that object's own reference-count
word**. That distinction decides the remedy, so the measurement is set out in
full below before the fix.

---

## 1. The measurement

### 1.1 A benchmark with nothing shared in it

`benchmarks/rcp_scaling.cpp` (new). N threads, each building and tearing down
expression trees over *its own* symbols:

~~~c++
RCP<const Basic> a = mul(x, y);      // Mul, coefficient = one
RCP<const Basic> b = add(a, z);      // Add, coefficient = zero
RCP<const Basic> c = sub(b, x);      // multiplies by minus_one
RCP<const Basic> d = pow(c, integer(2 + (i & 1)));
RCP<const Basic> e = div(d, y);
~~~

Two threads share no expression, no symbol and no dictionary. The only thing
they have in common is the library's process-global singletons — and nearly
every node built takes a reference to one of them, because `Mul` stores `one`
in its coefficient slot, `Add` stores `zero`, and `sub()` multiplies by
`minus_one`.

Release, `cooperative_intrusive` + `WITH_SYMENGINE_THREAD_SAFE`, 32-core box,
300 000 iterations/thread, best of 3:

| threads | seconds | Mops/s | speedup |
|--------:|--------:|-------:|--------:|
| 1 | 0.145 | 2.073 | 1.00× |
| 2 | 0.404 | 1.483 | 0.72× |
| 4 | 1.708 | 0.703 | 0.34× |
| 8 | 4.168 | 0.576 | 0.28× |
| 16 | 9.341 | 0.514 | 0.25× |

Aggregate throughput *falls* fourfold as cores are added. The corpus workload
that produced the 2.02× figure in the ledger is doing plenty of other work
alongside; this benchmark is nearly pure singleton traffic, so it shows the
same defect undiluted. (§3.1 repeats this sweep at a longer setting, on both
backends, against the fixed tree.)

### 1.2 `perf c2c` names the lines

`perf_event_paranoid` was `-1` and `perf c2c` worked. Recording the 8-thread
run:

~~~
Load Local HITM                   :         77
Load Remote HITM                  :         27
Total Shared Cache Lines          :          3      <- 104 of 104 HITMs
~~~

| line | share of HITMs |
|---|---|
| `0x…074300` | 87.50% |
| `0x…074380` | 8.65% |
| `0x…0742c0` | 3.85% |

The three addresses are 0x40 and 0xc0 apart. Printing the singleton addresses
from a separate probe program gives the same relative layout — `zero` on line
L at offset 0x10, `one` on L+0x40 at offset 0x30, `minus_one` on L+0xc0 at
offset 0x10 — so the three hot lines are exactly **`zero`, `one` and
`minus_one`**. (`two`, on L+0x100, never appears.)

Within each line the pareto concentrates on one offset, and it is the object's
reference count. `Basic` derives from `EnableRCPFromThis`, so the layout is
vptr at +0 and `refcount_` at +8:

| line | object at | hot offset | = |
|---|---|---|---|
| L+0x40 | 0x30 | **0x38** | `one->refcount_` |
| L+0xc0 | 0x10 | **0x18** | `minus_one->refcount_` |
| L | 0x10 | **0x18** | `zero->refcount_` |

Offset 0x38 alone carries ~85% of the HITMs on the busiest line, i.e. roughly
three quarters of the whole run's. The instruction addresses are
`symengine_rcp.h:90/96/111` — the counter's relaxed load and its two
compare-exchanges — reached from `mul`, `Mul::Mul`, `Mul::from_dict`,
`Add::dict_add_term`, `Add::as_coef_term`, `sub`, `pow`, and from `RCP`
destructors inside `std::_Hashtable` and `std::_Rb_tree` nodes. That is the
ledger's list, confirmed independently.

**So the mechanism is: every thread performs an atomic read-modify-write on
the same three words, millions of times a second, and the cache-coherence
protocol serialises them.** Each hot line has exactly one hot word.

---

## 2. Remedy

### 2.1 Chosen: immortal reference counts

An object that is never destroyed does not need its reference count kept — no
sequence of increments and decrements can change the answer to "is it still
alive?". So don't keep it. `mark_immortal()` puts the counter into a state
that the reference operations recognise and **decline to write**, and
`constants.cpp` marks every `DEFINE_CONSTANTS` singleton once they are all
built.

The line then stays `Shared` in every core's cache instead of being dragged
`Modified` into one of them, which eliminates the traffic rather than
relocating it.

Encoding, per backend. Immortality is a **band** of counts tested by range,
not a distinguished value tested by equality — §2.5 is the whole argument for
why, and it is an ABI argument, not an aesthetic one. `mark_immortal()` writes
the band's centre in both intrusive backends.

* **`cooperative_intrusive`** — the state word stays odd, so the object
  remains "C++-owned" for every existing predicate, and its count field
  (`v >> 1`) is put in the band `[2^61, 3*2^61)` with the mark at `2^62`.
  Spelled out in `symengine_rcp.h` as `cooperative_immortal_lo`,
  `cooperative_immortal` and `cooperative_immortal_hi` — `(1<<62)|1`,
  `(1<<63)|1` and `(3<<62)|1` for a 64-bit `uintptr_t` — with
  `cooperative_state_is_immortal()` the range test.
  `is_uniquely_owned_by_cpp()` reports false, which is the right answer to
  "may I steal this?" about a singleton.
* **`symengine` (default)** — the bare `unsigned int` count is put in the band
  `[2^30, 3*2^30)` with the mark at `2^31` (`detail::immortal_refcount_lo`,
  `immortal_refcount`, `immortal_refcount_hi` and `refcount_is_immortal()`).
  `RCP<T>` now goes through `inc_ref()` / `dec_ref()` on `EnableRCPFromThis`
  instead of touching `refcount_` directly, which is how the cooperative
  backend's `RCP<T>` was already written.
* **`teuchos`** — keeps its count in a separate node rather than in the
  object, has no equivalent, and `mark_immortal()` is a documented no-op
  there. Teuchos is not the thread-scalable backend in the first place.

**`UINT_MAX` survives in this design, but only as a report.** `use_count()`
maps any in-band value to `~0u` in both intrusive backends
(`symengine_rcp.h`'s `EnableRCPFromThis::use_count()` and
`symengine_cooperative_intrusive_counter::use_count()` in
`symengine_rcp_cooperative.cpp`), because no count is kept and "effectively
infinite" is the only honest answer. That is a **reporting convention**, and
it is deliberately not the storage encoding — an earlier version of this note
conflated the two, and the conflation is what made the retired top-of-range
sentinel look safe.

One observable consequence of keeping them separate: a consumer compiled
**before** the masking existed reads `refcount_` inline and so sees the raw
biased count through the public API. Measured against a current library:
`use_count=2147483648 is_immortal=0 value=1`. Harmless — that consumer's
reference traffic is RAII-balanced, so the count hovers at the band centre and
the object stays alive — but observable, and worth knowing before someone
reports it as corruption.

### 2.2 Rejected: cache-line padding of the singleton storage

Padding was the obvious cheap option and the ledger row suggests it
("cache-line padding of the global singletons is not [handled]"). The `perf
c2c` pareto is what rules it out: **each of the three lines has one hot word,
and it is that object's own count.** Padding separates `one` from whatever
shares its line — on line L+0x40 that is a GMP limb allocation, which does
contribute reads — but it cannot separate `one`'s count from itself. Sixteen
threads incrementing `one->refcount_` contend on that word whether or not it
has 56 bytes of padding around it. Best case padding converts one badly-shared
line into three separately-shared ones; the measured ceiling is set by the
count, not by its neighbours.

It is also more invasive than it sounds: the constants are ordinary heap
objects from `make_rcp`, so padding them means either over-aligning every
`Basic` (a large memory cost on every expression node in the system) or giving
the singletons a bespoke allocation path.

### 2.3 Costs of what was chosen

* **A branch on the hot path.** One load-and-compare before the atomic
  read-modify-write. Told nothing, GCC lays the immortal early-return out as
  the fallthrough and the atomic as a cold block, which cost a consistent
  2–3% single-threaded; a `__builtin_expect` that reference operations are
  usually on mortal objects puts the increment back inline and brings that
  down to the ~1% measured in §3.2. The macro applies the bias directly to the
  tested expression on both supported compiler families; it is not a
  language-floor compatibility claim.
* **Memory: the constants are never freed.** A few dozen objects plus the
  small expression trees below them (`sq3`, `C0`, `mC0`, …) — order 40 kB.
  They stay reachable from the static storage that names them, so a leak
  checker sees *still reachable*, not *leaked*; §4.3 records the ASan run.
* **No layout change, and a *checked* ABI claim.** No type changes size or
  alignment; `static_assert(sizeof(counter) == sizeof(void *))` still holds.
  The state word's encoding gains one tag bit (§2.4), which is internal to the
  counter class.

  That is not sufficient on its own, and the first version of this note
  claimed too much from it. Reference traffic is *inline*, so a consumer built
  against an older header runs its own arithmetic on objects this library
  marks — and the SONAME does not change, so replacing the shared object
  presents that pairing as compatible. With immortality first written as a
  top-of-range sentinel it was not: the old inline increment wrapped `~0u` to
  zero and the next destruction deleted a process-global singleton
  (heap-use-after-free under ASan; `after-a use_count=0`, `after-b
  use_count=1` without it). The sentinel is therefore a *band* of biased
  counts, chosen so that old inline arithmetic is merely un-optimised — see
  §2.5 — and `bin/test_abi_old_headers.sh` compiles a consumer against the
  pre-immortality headers, links it to the current library and runs it under
  ASan, on every CI run. No lane that builds library and consumer from one
  tree can see this class of break.
* **One tightened precondition**, fork-local: `set_self_external()` now
  requires a 4-byte-aligned foreign pointer where two bytes used to do. It is
  checked with a loud abort rather than assumed. See §4.2.

### 2.4 The part that had to be got right: the foreign hand-off

The fork's whole reason for existing is the cooperative-intrusive backend,
where the state word is a union: bit 0 set means a C++ reference count, bit 0
clear means a pointer to a foreign-runtime object that owns this one.
`set_self_external(o)` performs the transition by folding the current C++
count into the wrapper's own count — `incref(o)` once per outstanding C++
reference — so that every later C++ decrement is one of the wrapper's
decrements.

An immortal object has no count to fold. This is not a detail that can be
waved past: the bindings **do** wrap the singletons —
`nbsymengine/support/nanobind_module_common.h` lists `zero`, `one`,
`minus_one`, `two`, `I`, `pi`, `E`, the `Inf` family, `Nan`, the booleans and
the whole trig table as objects it calls `set_self_external()` on at module
import. A design in which immortality and externalisation are mutually
exclusive would break every cooperative binding, and one that transferred a
*guessed* count would under-reference the wrapper by however many `Mul`
coefficient slots happen to point at `one` at that moment, freeing it early.

So immortality survives externalisation. Bit 1 of the external state word
marks it:

| state word | meaning |
|---|---|
| bit 0 set, outside the band | C++-owned, count = `v >> 1` — unchanged |
| bit 0 set, inside the band | immortal, no wrapper |
| bit 0 clear, bit 1 clear | external-owned — unchanged |
| bit 0 clear, bit 1 set | immortal, wrapper attached |

("the band" is `[cooperative_immortal_lo, cooperative_immortal_hi)`; see §2.1
and §2.5. An immortal *external* word is recognised by bit 1 rather than by
range, because bits 2 and up are a pointer there, not a count.)

`set_self_external()` on an immortal object grants the wrapper **exactly one
permanent reference** and stays immortal. That one reference is what keeps
`self_external()` from ever returning a dangling pointer: C++ can reach the
object forever, so the wrapper must live forever too. C++ reference traffic
then stays off both the hooks and the state word — so a Python or Perl or PHP
process keeps the scaling as well, which the alternative designs would have
given up.

Everything a binding can observe behaves as it did:

* `self_external()` returns the wrapper (bit 1 masked off).
* `is_external_owned()` is true.
* `use_count()` is 0, as for any externalised object.
* `is_uniquely_owned_by_cpp()` is false, so no move-out path steals from it.
* `detach_external()` returns the wrapper and restores
  immortal-*without*-a-wrapper — not a count of zero, which would make the
  next `dec_ref()` an underflow abort.

`nbsymengine`'s shutdown protocol (`register_singleton_cleanup`) still works
unmodified: it detaches, reads `R = Py_REFCNT(o)`, calls `inc_ref()` `R+2`
times and `Py_DECREF`s `R` times. The `inc_ref()` calls become no-ops, the
`Py_DECREF`s release the wrapper including the one permanent reference, and
the C++ object stays immortal. Its careful accounting is *vacuous* under
immortality rather than wrong — it does not depend on how many references were
transferred at hand-off, only on the wrapper's Python count at shutdown.

### 2.5 The other part that had to be got right: the sentinel's value

Immortality is a *band* of reference-count values, not one distinguished
value, and the reason is the ABI point in §2.3: reference traffic is inline,
so the arithmetic that runs on a marked object is whatever header the consumer
was compiled with, and a consumer compiled before immortality existed does not
know the mark is there.

| count | meaning |
|---|---|
| `< lo` | ordinary, mortal, counted by everybody |
| `[lo, hi)` | immortal — this header elides the read-modify-write |
| `>= hi` | ordinary again |

`mark_immortal()` writes the band's centre. `lo`, the centre and `hi` are at
2^30, 2^31 and 3·2^30 for the original backend's 32-bit count, and at 2^61,
2^62 and 3·2^61 for the cooperative state word's count field.

An old consumer's traffic is RAII-balanced, so the count hovers at the centre
and never approaches zero: that consumer is correct, merely un-optimised,
which is what it had before immortality existed. A current consumer takes the
early return and elides the write, so §3.1's scaling is unaffected.

The objection is drift, because a *mixed* process can produce unbalanced
operations: a reference the library creates (elided) and the old consumer
destroys (counted) is a net decrement, and the mirror image a net increment.
Leaving the band takes 2^30 of them for the 32-bit count. But the size of that
margin is the lesser argument — **both edges are absorbing**. The moment a
count leaves the band, this header counts truthfully again, so the very
traffic that was unbalanced becomes balanced and the count cannot walk
further. Below the low edge it can fall only by the number of *simultaneously
live* references whose increment was elided, and 2^30 live references is 8 GB
of pointers alone.

The cost is one that the previous sentinel also had, moved down by a factor of
four: an ordinary object holding `lo` or more references would be misread as
immortal. 2^30 references cannot fit in an address space that must also hold
the objects making them.

The alternatives were to bump the SONAME, or to keep the mark off exported
constants. The first is correct and was the fallback; it was not chosen
because it forces every consumer of a fork whose singletons are the entire
point to rebuild, in exchange for a guarantee the band already provides. The
second gives up the optimisation on exactly the objects it exists for.

---

## 3. Before and after

### 3.1 Scaling

Same box, Release, 400 000 iterations/thread, best of 7, base and fixed trees
run back to back in the same session. "base" is commit `5f088693` (the
benchmark alone, before any counter change), "after" is `HEAD`. The box is
shared with another agent; the run below was taken at a 1-minute load average
of 2.3, and the whole sweep was repeated in three separate sessions at
different loads, with the 16-thread speedups landing at 0.24–0.25× / 4.4–4.6×
(cooperative) and 0.37–0.38× / 4.1–4.2× (default) every time.

**`cooperative_intrusive` + `WITH_SYMENGINE_THREAD_SAFE`** — the
configuration `symengine.rubi` uses:

| threads | before Mops/s | before | after Mops/s | after |
|--------:|--------------:|-------:|-------------:|------:|
| 1 | 2.185 | 1.00× | 2.135 | 1.00× |
| 2 | 1.625 | 0.74× | 4.061 | 1.90× |
| 4 | 1.213 | 0.55× | 7.512 | 3.52× |
| 8 | 0.787 | 0.36× | 10.077 | 4.72× |
| 16 | 0.518 | 0.24× | 9.884 | 4.63× |

**`symengine` (default) + `WITH_SYMENGINE_THREAD_SAFE`** — what upstream
would see:

| threads | before Mops/s | before | after Mops/s | after |
|--------:|--------------:|-------:|-------------:|------:|
| 1 | 2.560 | 1.00× | 2.360 | 1.00× |
| 2 | 2.143 | 0.84× | 4.439 | 1.88× |
| 4 | 1.733 | 0.68× | 7.291 | 3.09× |
| 8 | 0.862 | 0.34× | 7.716 | 3.27× |
| 16 | 0.982 | 0.38× | 9.705 | 4.11× |

At 16 threads that is **19.1×** more aggregate throughput on the cooperative
backend and **9.9×** on the default one.

`perf c2c` on the fixed build, same 8-thread workload:

~~~
Total Shared Cache Lines          :          1
Load Local HITM                   :          1
~~~

Down from 3 lines and 104 HITMs to one line and one HITM — and that one is a
*read* of the memoised `hash_` word, not a count.

The remaining sublinearity above 4 threads is the allocator, not coherence: a
16-thread CPU profile of the fixed build is `mul` 12.5%, hashtable assign
11.2%, `Mul::Mul` 9.3%, `Mul::~Mul` 5.9% — real work — with `cfree` 6.4%,
`_int_malloc` 4.5% and `malloc` 3.5% behind them. That is a separate problem
and not this defect.

### 3.2 Single-threaded cost

Wall clock could not resolve this on a shared box — repeated pinned runs
disagreed by more than the effect. Retired instructions can, being exact and
load-independent, so that is what is reported; cycles are given alongside for
the one case where the box was quiet enough to mean anything.

| workload | backend | instructions before → after | |
|---|---|---|---|
| `expand1` | symengine | 541.9 M → 550.6 M | **+1.6%** |
| `expand1` | cooperative | 564.0 M → 572.9 M | **+1.6%** |
| `rcp_scaling`, 1 thread | symengine | 4 254.9 M → 4 479.7 M | +5.3% |
| `rcp_scaling`, 1 thread | cooperative | 4 926.4 M → 4 983.8 M | +1.2% |

`rcp_scaling` is deliberately singleton-saturated, so its +5.3% is an upper
bound rather than a representative figure; `expand1` is the realistic one.

Cycles tell a better story than instructions do, because the check *replaces*
a locked read-modify-write whenever it fires, and a `lock addl` costs far more
than a load and a compare. Four alternating pinned `expand1` runs, default
backend:

| | base | after |
|---|---|---|
| cycles | 308.7 / 302.3 / 292.6 / 287.6 M | 297.8 / 297.1 / 303.0 / 286.5 M |
| mean | 297.8 M | 296.1 M |

i.e. **cycle-neutral within noise** on a realistic single-threaded workload,
at +1.6% instructions. Before the branch-layout hint of commit `e00dcc88` it
was a consistent 2–3% slower in wall clock.

---

## 4. Correctness

### 4.1 Both backends, and the thread-safety flag

Every configuration below is Release unless noted, `BUILD_TESTS=ON`,
`ctest --output-on-failure -E "test_parser"`:

| backend | `WITH_SYMENGINE_THREAD_SAFE` | result |
|---|---|---|
| `cooperative_intrusive` | ON | 61/61 |
| `cooperative_intrusive` | OFF | 61/61 |
| `symengine` | ON | 61/61 |
| `symengine` | OFF | 61/61 |
| `teuchos` | ON | 61/61 |
| `cooperative_intrusive`, `-std=c++11` | ON | 61/61 |
| `cooperative_intrusive`, `-std=c++14` | ON | 61/61 |
| `symengine`, `-std=c++11` | ON | 61/61 |

(The sanitizer builds below register 62 rather than 61 because they use
`INTEGER_CLASS=boostmp`, which adds `test_integer_class`.)

**On the language standard.** The three `-std=` rows above are real: this
branch declares `target_compile_features(symengine PUBLIC cxx_std_11)`, so
the pre-C++17 sides of the two `__cplusplus` tests in `symengine_rcp.h` are
reachable configurations and were built and run, not merely parsed. Nothing in
the immortality change needs a newer language: `immortal_refcount` is an
enumerator and cannot require storage, and the band predicates are ordinary
unsigned arithmetic.

Two pre-existing C++-level notes, neither introduced by this series and both
visible only below the compiler default of C++17:

* `cooperative_intrusive_init()`'s hook types carry `noexcept` only from
  C++17, and since C++17 that is part of the type, so the function's mangled
  name differs between a C++11/14 and a C++17 consumer (`PFvPvE` against
  `PDoFvPvE`). A binding must be compiled on the same side of that test as
  the library. This is now stated at the declaration.
* `test_cooperative_intrusive_rcp.cpp`'s race test uses a lambda init-capture
  to move an `RCP` into the worker thread, which is C++14. GCC and Clang
  accept it under `-std=c++11` with a warning; a strict C++11 compiler would
  not.

New tests:

* `test_rcp.cpp` — `mark_immortal()` per backend: the count stops moving,
  the object survives its last `RCP`, `is_uniquely_owned()` goes false, and
  the library's own constants come back `is_immortal()`.
* `test_cooperative_intrusive_rcp.cpp` — an immortal object keeps no count and
  reports `UINT_MAX`; the hand-off grants exactly one foreign reference (not
  one per outstanding reference), keeps later C++ traffic off the hooks
  entirely, and `detach_external()` round-trips to immortal-without-a-wrapper
  so the final `RCP` release is not an underflow.

`mark_immortal()` itself is a compare-exchange on the cooperative backend, so
it composes with concurrent reference operations: a thread whose exchange
loses to it retries, sees the immortal state and stops, silently abandoning a
reference that no longer needs counting. What it cannot survive is racing the
*last* `dec_ref()`, because that thread has already been told to delete; so
the caller must hold a reference, which `constants.cpp` does (the global `RCP`
that names the constant) and which the doc comment states.

**Two transitions were not linearizable, and both are now.** On the original
intrusive backend `inc_ref()` was a test and a *separate* read-modify-write,
with marking a plain store between them; the increment then moved the count
off the mark. Each operation is individually atomic, so this is a lost update
rather than a data race and **TSan cannot flag it** — and the tests above
generate traffic only *after* marking, never across the transition. Both
`inc_ref()` and `dec_ref()` now re-test the value the read-modify-write
returns and undo it if marking won; a compare-exchange loop would also work
but would cost every ordinary object a retry loop to linearize a transition
that happens a few dozen times per process, at static init.

On the cooperative backend `set_self_external()` recognised marking on its
initial load but not if marking won between that load and the publishing
exchange: the retry read the sentinel as an ordinary count and replayed it
into the foreign runtime, ~2^62 increfs, i.e. a hang. Immortality is now
re-tested after *every* failed exchange, and recognising it switches to the
one-permanent-reference protocol rather than resuming the count-folding one.
The regression for it is deterministic, not a stress test: the foreign incref
hook blocks on its first call, which parks externalization exactly in the
window, the mark is made from the test thread, and only then is the hook
released. The hook aborts past a small bound, so the pre-fix behaviour is a
failure in milliseconds instead of a test run that never ends.

### 4.2 The tightened precondition

`set_self_external()` now requires `((uintptr_t)o & 3) == 0`, where before
only bit 0 mattered. A pointer with bit 1 set would read back as immortal and
every reference operation on that object would be silently dropped — a
use-after-free — so it is checked with an abort rather than assumed. Every
runtime the backend targets satisfies it: `PyObject` (16-byte aligned),
`zend_object`, Perl `SV`, Swift class instances are all at least
pointer-aligned. The one place that did *not* satisfy it was the existing test
fixture, which used `alignas(2)` to document the old contract; it is now
`alignas(4)`.

This is a fork-local API. Upstream has no `set_self_external` and this
paragraph does not apply to the upstream PR.

### 4.3 Sanitizers

**TSan** (clang, instrumented libc++, `cooperative_intrusive` +
thread-safe, Debug): 61/62 pass. The one failure is `test_functions`, which
reports six data races, all of them inside `std::cout`'s `ios_base` flags and
width, written concurrently by the sixteen threads of `test_dummy` from a
leftover `std::cout << "out=" << out` debug line in the test's own thread body
(`test_functions.cpp:4218`). No frame in any of the six stacks is in
reference-counting code. Reproduced identically on base commit `5f088693`; the
file is not touched by this branch. Pre-existing test defect, worth a separate
one-line fix.

The benchmark itself was also run under TSan at 2, 4 and 8 threads — a
workload whose entire point is hammering the immortal counters from several
threads at once — with no warnings.

The immortality check does not touch the memory ordering. In the cooperative
backend it is a test on the state word placed *before* the release
compare-exchange and the acquire fence, both untouched; and the value it tests
cannot transition under it, because `mark_immortal()` is one-way. In the
`symengine` backend the added load is sequentially consistent, like the
increment it guards, so nothing was weakened.

**ASan + LeakSanitizer** (clang, instrumented libc++, `cooperative_intrusive`
+ thread-safe, Debug, `detect_leaks=1`, no suppressions): **62/62**.

That is the number that matters for §2.3's memory claim: the library's
constants are now retained for the life of the process, and LSan does not
report them, because they stay reachable from the static storage that names
them. The first ASan run did fail, on the new `test_rcp` case — which had
marked an object immortal and then dropped the only pointer to it off a stack
frame. That is a real leak and LSan was right; commit `52a4f8b6` parks the
pointer in a static, which is the actual model of what the constants are.

### 4.4 What could not be checked here

* The foreign bindings themselves (`nbsymengine`, `symengine.pl`,
  `symengine.php`, `symengine.swift`) live outside this worktree and were not
  rebuilt against this library. §2.4 is an argument from their source, not a
  test run. Root CI lanes ci-02 … ci-08 are the confirmation, and ci-04 (the
  leak lane) is the one to watch, because the constants are now retained.
* The `symengine.rubi` corpus is the supervisor's separate confirmation, per
  the brief.

---

## 5. For the upstream PR

The generic core is three files and no new build option:

* `symengine/symengine_rcp.h` — `mark_immortal()` / `is_immortal()` on
  `EnableRCPFromThis`, the immortal *band* per backend (§2.5 — not a
  sentinel, and the difference is the ABI argument), the branch in the
  reference operations, and `RCP<T>` (symengine backend) routed through
  `inc_ref()`/`dec_ref()` instead of touching `refcount_`.
* `symengine/constants.cpp` — mark the `DEFINE_CONSTANTS` singletons, after
  they are all built.
* `benchmarks/rcp_scaling.cpp` + its CMake entry — the evidence, reproducible.
* tests in `symengine/tests/rcp/test_rcp.cpp`.

Everything about the cooperative-intrusive backend — the bit-1 external tag,
the one-permanent-reference hand-off, the alignment precondition — is
fork-local and drops out of the upstream diff entirely, because upstream has
no such backend.

Three points a reviewer will reasonably raise:

1. **"The constants leak now."** They are retained, and reachable from the
   statics that name them, so leak checkers class them as still-reachable —
   the full suite is clean under LSan with `detect_leaks=1` and no
   suppressions. That is the trade: a few dozen objects held to the end of the
   process, in exchange for the reference counts of the most-referenced
   objects in the library never being written. The same trade CPython made in
   PEP 683.
2. **"A branch on the hottest path in the library."** +1.6% instructions and
   no measurable cycle cost on `expand1`, against 10–19× multi-threaded,
   because the check replaces a locked read-modify-write whenever it fires.
   If a reviewer wants it gated, the natural gate is a CMake option that makes
   `mark_immortal()` a no-op everywhere — but these numbers do not seem to
   justify the configuration surface.
3. **"Why not just pad?"** §2.2, with the `perf c2c` pareto: each hot line
   has one hot word and it is that object's own count.

Two adjacent cleanups this change makes visible but deliberately does *not*
do, to keep the diff reviewable:

* After routing `RCP<T>` through `inc_ref()`/`dec_ref()`, the `symengine` and
  `cooperative_intrusive` specialisations of `RCP<T>` are textually
  identical. They could be merged, deleting about a hundred duplicated lines.
* `Basic::hash()` memoises through an explicitly relaxed atomic; the same
  argument applies to the reference count's increment, which is currently
  sequentially consistent in the `symengine` backend where relaxed would do
  (the decrement genuinely needs release, as the cooperative backend's comment
  already explains). That is an ordering change and belongs in its own PR.

---

## 6. Commits

| | |
|---|---|
| `5f088693` | `benchmarks: add rcp_scaling` — the evidence, before any fix |
| `a85f6b0b` | `rcp: give the process-global singletons immortal reference counts` |
| `e00dcc88` | `rcp: hint that reference operations are usually on mortal objects` |
| `b5041050` | `rcp: check the foreign pointer's alignment instead of assuming it` |
| `52a4f8b6` | `tests: park the immortal test objects in statics` |
| `b997aa79` | `docs:` this file |
| `957de7ea` | `rcp: shorten the immortal-state names and restore clang-format` |
| `918a04fa` | `rcp: spell out use_count() for an immortal object` |
| *(and the commit that adds this table)* | |

These are this branch's revisions. The series was developed on the fork's
`rubi` branch and cherry-picked here; each commit carries a
`(cherry picked from commit ...)` trailer naming its original, and the
follow-ups that came after this table was first written -- the immortal band,
the two linearizability fixes, the acquire on the final decrement, and
`bin/test_abi_old_headers.sh` -- are in the log above this one.

`a85f6b0b` is the change; the ones after it are things found while checking
it, each of which would have been a defect had it shipped.

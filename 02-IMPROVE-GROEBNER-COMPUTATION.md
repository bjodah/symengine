# Background

We have (hopefully) implemented 01-IMPL-PLAN-GROEBNER.md at this point.

# Goals

1. I would like to add some relevant benchmarks (to the ./benchmarks/ folder). I've added one suggested idea here:
```console
python demo_groebner_bench_sympy.py 
timing1: 0.7096274869982153
timing2a: 0.014813994988799095
timing2b: 1.2804055586457253e-05
```
If you can construct more benchmarks that would be wonderful, SymPy's source code (which includes distributed tests) might be a good source (search for fglm, groebner and f5b).

2. Currently, we also have MoGVW method listed, but not implemented, I've added the LaTeX source of the arxiv paper as `./mogvw.tex`. Another method I'd like us to add is "M4GB", I've added the reference implementation git repo under `./m4gb`. I want you to write a highly detailed implementation plans for implementing both MoGVW and M4GB. Save the plan as `03-IMPL-PLAN-MOGVW-AND-M4GB.md`.

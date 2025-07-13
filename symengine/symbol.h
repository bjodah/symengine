/**
 *  \file symbol.h
 *  Class Symbol
 *
 **/
#ifndef SYMENGINE_SYMBOL_H
#define SYMENGINE_SYMBOL_H

#include <symengine/basic.h>

namespace SymEngine
{

class Symbol : public Basic
{
private:
    //! name of Symbol
    std::string name_;

public:
    IMPLEMENT_TYPEID(SYMENGINE_SYMBOL)
    //! Symbol Constructor
    explicit Symbol(const std::string &name);
    //! \return Size of the hash
    hash_t __hash__() const override;
    /*! Equality comparator
     * \param o - Object to be compared with
     * \return whether the 2 objects are equal
     * */
    bool __eq__(const Basic &o) const override;

    /*! Comparison operator
     * \param o - Object to be compared with
     * \return `0` if equal, `-1` , `1` according to string compare
     * */
    int compare(const Basic &o) const override;
    //! \return name of the Symbol.
    inline const std::string &get_name() const
    {
        return name_;
    }

    vec_basic get_args() const override
    {
        return {};
    }
    RCP<const Symbol> as_dummy() const;
};

class Dummy : public Symbol
{
private:
    //! Dummy count
#ifdef WITH_SYMENGINE_THREAD_SAFE
    static size_t count_;
    // static std::atomic<size_t> count_;
#else
    static size_t count_;
#endif
    //! Dummy index
    size_t dummy_index;

public:
    IMPLEMENT_TYPEID(SYMENGINE_DUMMY)
    //! Dummy Constructors
    explicit Dummy();
    explicit Dummy(const std::string &name);
    Dummy(const std::string &name, size_t index);
    //! \return Size of the hash
    hash_t __hash__() const override;
    /*! Equality comparator
     * \param o - Object to be compared with
     * \return whether the 2 objects are equal
     * */
    bool __eq__(const Basic &o) const override;
    /*! Comparison operator
     * \param o - Object to be compared with
     * \return `0` if equal, `-1` , `1` according to string compare
     * */
    int compare(const Basic &o) const override;
    size_t get_index() const
    {
        return dummy_index;
    }
    static constexpr const char *default_Dummy_prefix_ = "_Dummy_";
    static constexpr size_t default_Dummy_prefix_len_
        = sizeof(default_Dummy_prefix_) - 1;
    static_assert(default_Dummy_prefix_len_ == 7);
};

//! inline version to return `Symbol`
inline RCP<const Symbol> symbol(const std::string &name)
{
    return make_rcp<const Symbol>(name);
}

//! inline version to return `Dummy`
inline RCP<const Dummy> dummy()
{
#ifdef WITH_SYMENGINE_THREAD_SAFE
    return make_rcp<const Dummy>();
    // return make_rcp<const Dummy>(Dummy::default_Dummy_prefix_,
    // Dummy::count_.fetch_add(1, std::memory_order_relaxed));
#else
    return make_rcp<const Dummy>();
#endif
}

inline RCP<const Dummy> dummy(const std::string &name)
{
    return make_rcp<const Dummy>(name);
}

inline RCP<const Dummy> dummy(const std::string &name, size_t dummy_index)
{
    return make_rcp<const Dummy>(name, dummy_index);
}

inline bool is_a_Symbol(const Basic &b)
{
    return b.get_type_code() == SYMENGINE_SYMBOL
           || b.get_type_code() == SYMENGINE_DUMMY;
}

} // namespace SymEngine

#endif

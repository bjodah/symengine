#ifndef SYMENGINE_PARSER_H
#define SYMENGINE_PARSER_H

#include <symengine/basic.h>

#include <memory> // std::shared_ptr

namespace SymEngine
{

struct ParserSettings {
    static std::shared_ptr<ParserSettings> make_default();

    bool convert_xor{true};
    std::shared_ptr<std::map<std::string, RCP<const Basic>>> constants{};
    std::shared_ptr<
        std::map<std::string, std::function<RCP<const Basic>(vec_basic &)>>>
        multi_arg_functions{};
    std::shared_ptr<std::map<
        std::string, std::function<RCP<const Basic>(const RCP<const Basic> &,
                                                    const RCP<const Basic> &)>>>
        double_arg_functions{};
    std::shared_ptr<std::map<
        std::string, std::function<RCP<const Basic>(const RCP<const Basic> &)>>>
        single_arg_functions{};
    std::shared_ptr<std::map<
        std::string, std::function<RCP<const Basic>(const RCP<const Basic> &)>>>
        tripple_arg_functions{};
};

std::shared_ptr<std::map<std::string, RCP<const Basic>>>
get_default_parser_constants();
std::shared_ptr<std::map<
    std::string, std::function<RCP<const Basic>(const RCP<const Basic> &,
                                                const RCP<const Basic> &)>>>
get_default_double_arg_functions();
std::shared_ptr<
    std::map<std::string, std::function<RCP<const Basic>(vec_basic &)>>>
get_default_multi_arg_functions();
std::shared_ptr<std::map<
    std::string, std::function<RCP<const Basic>(const RCP<const Basic> &)>>>
get_default_single_arg_functions();

RCP<const Basic> parse(const std::string &s,
                       std::shared_ptr<const ParserSettings> settings = {});
RCP<const Basic> parse(const std::string &s, bool convert_xor);
RCP<const Basic>
parse(const std::string &s, bool convert_xor,
      const std::map<const std::string, const RCP<const Basic>> &constants);

RCP<const Basic> parse_old(const std::string &s, bool convert_xor = true);
RCP<const Basic>
parse_sbml(const std::string &s,
           const std::map<const std::string, const RCP<const Basic>> &constants
           = {});
} // namespace SymEngine

#endif

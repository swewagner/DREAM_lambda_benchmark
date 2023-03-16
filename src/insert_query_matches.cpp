#include <random>
#include <iostream>
#include <fstream>

#include <seqan3/argument_parser/all.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/io/sequence_file/output.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/utility/views/zip.hpp>

// struct for command line arguments
struct cmd_arguments
{
    std::filesystem::path query_file_path{};
    std::filesystem::path match_file_path{};
    std::filesystem::path out_file_path{};
};

struct dna4_input_trait : seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = seqan3::dna4;
};

void run_program(cmd_arguments const & args)
{
    // random number generator setup
    std::random_device rd;
    std::mt19937 gen(rd());

    // IO
    seqan3::sequence_file_input<dna4_input_trait> q_fin{args.query_file_path};
    seqan3::sequence_file_input<dna4_input_trait> m_fin{args.match_file_path};
    seqan3::sequence_file_output fout{args.out_file_path};

    using types = seqan3::type_list<std::vector<seqan3::dna4>, std::string>;
    using fields = seqan3::fields<seqan3::field::seq, seqan3::field::id>;
    using sequence_record_type = seqan3::sequence_record<types, fields>;

    for (auto && [q_rec, m_rec] : seqan3::views::zip(q_fin, m_fin))
    {
        // find pos to insert match
        std::uniform_int_distribution<uint64_t> match_start_dis(0, std::ranges::size(q_rec.sequence()) - std::ranges::size(m_rec.sequence()));
        uint64_t const match_start_pos = match_start_dis(gen);

        // insert match
        std::vector<seqan3::dna4> query_pre = q_rec.sequence() |
                                              seqan3::views::slice(0, match_start_pos) | seqan3::ranges::to<std::vector>();
        std::vector<seqan3::dna4> query_match = m_rec.sequence() |
                                                seqan3::ranges::to<std::vector>();
        std::vector<seqan3::dna4> query_suf = q_rec.sequence() |
                                              seqan3::views::slice(match_start_pos + std::ranges::size(m_rec.sequence()), std::ranges::size(q_rec.sequence())) | seqan3::ranges::to<std::vector>();
        std::vector<std::vector<seqan3::dna4>> v{ query_pre, query_match, query_suf };
        std::vector<seqan3::dna4> new_query = std::ranges::join_view(v) |
                                              seqan3::ranges::to<std::vector>();

        // write new query to output file
        sequence_record_type record{std::move(new_query), std::move(m_rec.id())};
        fout.push_back(record);
    }

}

void initialize_argument_parser(seqan3::argument_parser & parser, cmd_arguments & args)
{
    // add meta data
    parser.info.app_name = "insert_query_matches";
    // parser.info.author = "Swenja Wagner";
    // parser.info.email = "swenja.wagner@fu-berlin.de";
    // parser.info.date = "28.02.2023";
    parser.info.short_description = "inserts matches into given queries";

    // add options/flags
    parser.add_option(args.query_file_path,
                      'q',
                      "in_queries",
                      "Please provide a file with query sequences.",
                      seqan3::option_spec::required,
                      seqan3::input_file_validator{{"fa", "fasta", "fna", "ffn", "ffa", "frn"}});

    parser.add_option(args.match_file_path,
                      'm',
                      "in_matches",
                      "Please provide a file with match sequences.",
                      seqan3::option_spec::required,
                      seqan3::input_file_validator{{"fa", "fasta", "fna", "ffn", "ffa", "frn"}});

    parser.add_option(args.out_file_path, // TODO check possible outputs
                      'o',
                      "out",
                      "Please provide a file in which the query sequences with inserted matches should be written.",
                      seqan3::option_spec::required,
                      seqan3::output_file_validator{seqan3::output_file_open_options::open_or_create, {"fa", "fasta", "fna", "ffn", "ffa", "frn"}});
}

int main(int argc, char ** argv)
{
    // initialize parser
    seqan3::argument_parser command_line_parser{"match_insertion_parsing", argc, argv};
    cmd_arguments args{};

    initialize_argument_parser(command_line_parser, args);
    // add info, options, flags, pos options

    try
    {
        command_line_parser.parse(); // trigger command line parsing
    }
    catch(seqan3::argument_parser_error const & ext) // catch user errors
    {
        seqan3::debug_stream << "[Error] " << ext.what() << "\n";
        return -1;
    }

    run_program(args);

    return 0;
}
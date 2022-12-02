#include <random>

#include <seqan3/argument_parser/all.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/io/sequence_file/output.hpp>
#include <seqan3/core/debug_stream.hpp>

// struct for command line arguments
struct cmd_arguments
{
    std::filesystem::path in_file_path{};
    std::filesystem::path out_file_path{};
    uint32_t query_num{}; // TODO find suitable defaults
    uint32_t query_len{};
};

void initialize_argument_parser(seqan3::argument_parser & parser, cmd_arguments & args)
{
    // add meta data
    parser.info.app_name = "generate_query_matches";
    parser.info.author = "Swenja Wagner";
    parser.info.email = "swenja.wagner@fu-berlin.de";
    parser.info.date = "02.12.2022";
    parser.info.short_description = "generating query matches from given references";
    parser.info.version = "0.0.1";
    // parser.info.examples = {}

    // add options/flags
    parser.add_option(args.in_file_path, // TODO check possible inputs
                      'i',
                      "in",
                      "Please provide a file with reference sequences.",
                      seqan3::option_spec::required,
                      seqan3::input_file_validator{{"fa", "fasta"}});

    parser.add_option(args.out_file_path, // TODO check possible outputs
                      'o',
                      "out",
                      "Please provide a file in which the query sequences should be written.",
                      seqan3::option_spec::required,
                      seqan3::output_file_validator{seqan3::output_file_open_options::open_or_create, {"fa", "fasta"}});

    parser.add_option(args.query_num,
                      'n',
                      "num_of_queries",
                      "Please provide the number of query sequences, that should be generated");

    parser.add_option(args.query_len,
                      'l',
                      "len_of_queries",
                      "Please provide the length that the generated query sequences should have.");
}

void run_program(cmd_arguments const & args)
{
    // read reference file

    // distribute number of queries to reference seqs

    // go through references take random start positions for matches

    // write matches to output file

    // write file to keep track of reference match pairs 

}

int main(int argc, char ** argv)
{
    // initialize parser
    seqan3::argument_parser command_line_parser{"query_generation_parsing", argc, argv};
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
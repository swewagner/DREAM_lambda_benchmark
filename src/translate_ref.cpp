#include <seqan3/alphabet/views/translate.hpp>
#include <seqan3/argument_parser/all.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/io/sequence_file/output.hpp>
#include <seqan3/core/debug_stream.hpp>

// struct for command line arguments
struct cmd_arguments
{
    std::filesystem::path in_file_path{};
    std::filesystem::path out_file_path{};
};

struct dna4_input_trait : seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = seqan3::dna4;
};

void run_program(cmd_arguments const & args)
{
    // IO
    seqan3::sequence_file_input<dna4_input_trait> fin{args.in_file_path};
    seqan3::sequence_file_output fout{args.out_file_path};

    for (auto const & [sequence, id, quality] : fin) {
        
        std::vector<seqan3::aa27> prot_seq = sequence |
                        seqan3::views::translate_single |
                        seqan3::ranges::to<seqan3::aa27_vector>();

        // write sequence to output
        std::vector<seqan3::phred42> const prot_quality(prot_seq.size(), seqan3::assign_rank_to(40u, seqan3::phred42{}));
        fout.emplace_back(prot_seq, id, prot_quality);
    }
}

void initialize_argument_parser(seqan3::argument_parser & parser, cmd_arguments & args)
{
    // add meta data
    parser.info.app_name = "translate_ref_seqs";
    parser.info.author = "Swenja Wagner";
    parser.info.email = "swenja.wagner@fu-berlin.de";
    parser.info.date = "18.01.2022";
    parser.info.short_description = "translating reference sequences from dna to prot alphabet";
    parser.info.version = "0.0.1";
    // parser.info.examples = {}

    // add options/flags
    parser.add_option(args.in_file_path,
                      'i',
                      "in",
                      "Please provide a file with reference sequences in a dna alphabet.",
                      seqan3::option_spec::required,
                      seqan3::input_file_validator{{"fa", "fasta", "fna", "ffn", "ffa", "frn"}});

    parser.add_option(args.out_file_path,
                      'o',
                      "out",
                      "Please provide a file in which the translated reference sequences should be written.",
                      seqan3::option_spec::required,
                      seqan3::output_file_validator{seqan3::output_file_open_options::open_or_create, {"fa", "fasta", "fna", "ffn", "ffa", "frn"}});
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
#include <random>
#include <iostream>
#include <fstream>

#include <seqan3/argument_parser/all.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/io/sequence_file/output.hpp>
#include <seqan3/core/debug_stream.hpp>

// struct for command line arguments
struct cmd_arguments
{
    std::filesystem::path in_file_path{};
    std::filesystem::path out_file_path{};
    std::filesystem::path ground_truth_file_path{};
    uint32_t query_num{}; // TODO find suitable defaults
    uint32_t query_len{};
    uint32_t ref_num{};
    double max_error_rate{0.01};
    bool verbose_ids{false};
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
    seqan3::sequence_file_input<dna4_input_trait> fin{args.in_file_path};
    seqan3::sequence_file_output fout{args.out_file_path};

    // go through references take random start positions for matches
    uint32_t pos{0};
    std::ofstream ground_truth_file(args.ground_truth_file_path, std::ios_base::trunc);
    for (auto const & [sequence, id, quality] : fin)
    {
        pos ++;
        uint32_t num_of_matches_per_record = pos <= args.query_num % args.ref_num ? args.query_num / args.ref_num + 1 : args.query_num / args.ref_num;

        uint64_t const reference_length = std::ranges::size(sequence);
        std::uniform_int_distribution<uint64_t> match_start_dis(0, reference_length - args.query_len);
        std::vector<std::string> matches_per_record{};

        for (uint32_t current_match_number = 0; current_match_number < num_of_matches_per_record; ++current_match_number)
        {
            // extract match
            uint64_t const match_start_pos = match_start_dis(gen);
            std::vector<seqan3::dna4> match = sequence |
                                              seqan3::views::slice(match_start_pos, match_start_pos + args.query_len) |
                                              seqan3::ranges::to<std::vector>();

            // introduce errors
            uint32_t max_errors = (int) (args.query_len * args.max_error_rate);
            std::uniform_int_distribution<uint32_t> match_error_pos_dis(0, args.query_len -1);
            std::uniform_int_distribution<uint8_t> dna4_rank_dis(0, 3);
            for (uint8_t error_count = 0; error_count < max_errors; ++error_count)
            {
                uint32_t const error_pos = match_error_pos_dis(gen);
                seqan3::dna4 const current_base = match[error_pos];
                seqan3::dna4 new_base = current_base;
                while (new_base == current_base)
                    seqan3::assign_rank_to(dna4_rank_dis(gen), new_base);
                match[error_pos] = new_base;
            }

            // write matches to output file
            std::vector<seqan3::phred42> const quality(args.query_len, seqan3::assign_rank_to(40u, seqan3::phred42{}));
            std::string match_id = id + "_" + std::to_string(current_match_number);
            std::string meta_info{};

            if (args.verbose_ids)
            {
                meta_info += ' ';
                meta_info += "reference_file='" + std::string{args.in_file_path} + "'";
                meta_info += ", reference_id='" + id + "'";
                meta_info += ", start_position=" + std::to_string(match_start_pos);
                meta_info += ", errors=" + std::to_string(max_errors);
            }
            fout.emplace_back(match, match_id + meta_info, quality);
            matches_per_record.push_back(match_id);
        }

        // write matches to corresponding reference ID in ground truth file
        ground_truth_file << id; // reference id
        for (auto & elem : matches_per_record)
        {
            ground_truth_file << "," + elem; // match ids
        }
        ground_truth_file << "\n";
    }
    ground_truth_file.close();
// TODO: maybe it should be possible to have matches mapping to multiple references
// for that a script that inserts matches in references should be written
// the groundtruth needs to be adjusted accordingly
// TODO2: maybe randomly skip some references, so that some references do not have matches
}

void initialize_argument_parser(seqan3::argument_parser & parser, cmd_arguments & args)
{
    // add meta data
    parser.info.app_name = "generate_query_matches";
    parser.info.author = "Swenja Wagner";
    parser.info.email = "swenja.wagner@fu-berlin.de";
    parser.info.date = "28.02.2023";
    parser.info.short_description = "generating query matches from given references";
    parser.info.version = "0.0.1";
    // parser.info.examples = {}

    // add options/flags
    parser.add_option(args.in_file_path,
                      'i',
                      "in",
                      "Please provide a file with reference sequences.",
                      seqan3::option_spec::required,
                      seqan3::input_file_validator{{"fa", "fasta", "fna", "ffn", "ffa", "frn"}});

    parser.add_option(args.out_file_path, // TODO check possible outputs
                      'o',
                      "out",
                      "Please provide a file in which the query sequences should be written.",
                      seqan3::option_spec::required,
                      seqan3::output_file_validator{seqan3::output_file_open_options::open_or_create, {"fa", "fasta", "fna", "ffn", "ffa", "frn"}});

    parser.add_option(args.ground_truth_file_path,
                      'g',
                      "ground_truth",
                      "Please provide a file in which the ground truth should be written (txt).",
                      seqan3::option_spec::required);

    parser.add_option(args.query_num,
                      'n',
                      "num_of_queries",
                      "Please provide the number of query sequences, that should be generated");

    parser.add_option(args.query_len,
                      'l',
                      "len_of_queries",
                      "Please provide the length that the generated query sequences should have.");

    parser.add_option(args.ref_num,
                      'r',
                      "num_of_references",
                      "Please provide the number of reference sequences in the input file");
    
    parser.add_option(args.max_error_rate,
                      'e',
                      "max_error_rate",
                      "Please provide the max error rate for a query");

    parser.add_flag(args.verbose_ids,
                     'v',
                     "verbose-ids",
                     "Puts sampling information into the ID");
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
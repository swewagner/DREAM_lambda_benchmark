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
    std::filesystem::path query_file_path{};
    std::filesystem::path ref_file_path{};
    std::filesystem::path ref_out_dir{};
    std::filesystem::path ground_truth_file_path{};
    uint32_t match_len{};
    uint32_t ref_num{};
    uint32_t bin_num{};
    double max_error_rate{0.01};
    bool verbose_ids{false};
};

struct aa_input_trait : seqan3::sequence_file_input_default_traits_aa
{
    using sequence_alphabet = seqan3::aa27;
};

struct match_info
{
    std::string bin_id;
    std::string ref_id;
    std::string query_id;
    uint32_t match_len;
    uint64_t start_pos_in_ref;
    uint64_t start_pos_in_query;
};


using types = seqan3::type_list<std::vector<seqan3::aa27>, std::string>;
using fields = seqan3::fields<seqan3::field::seq, seqan3::field::id>;
using sequence_record_type = seqan3::sequence_record<types, fields>;

match_info insert_match(std::vector<std::vector<sequence_record_type>> & bins, std::vector<seqan3::aa27> match, uint32_t & bin_id, uint32_t seq_idx, u_int32_t match_len, std::string query_id, uint64_t query_match_start_pos)
{
    // random number generator setup
    std::random_device rd;
    std::mt19937 gen(rd());

    std::uniform_int_distribution<uint64_t> insert_start_dis(0, std::ranges::size(bins[bin_id][0].sequence()) - match_len);
    uint64_t insert_start_pos = insert_start_dis(gen);
    std::vector<seqan3::aa27> seq_pre = bins[bin_id][seq_idx].sequence() |
                                        seqan3::views::slice(0, insert_start_pos) | seqan3::ranges::to<std::vector>();
    std::vector<seqan3::aa27> seq_suf = bins[bin_id][seq_idx].sequence() |
                                        seqan3::views::slice(insert_start_pos + match_len, std::ranges::size(bins[bin_id][0].sequence())) | seqan3::ranges::to<std::vector>();
    std::vector<std::vector<seqan3::aa27>> v{ seq_pre, match, seq_suf };
    std::vector<seqan3::aa27> new_seq = std::ranges::join_view(v) |
                                        seqan3::ranges::to<std::vector>();

    // write new query to output file
    sequence_record_type new_record{std::move(new_seq), std::move(bins[bin_id][0].id())};
    bins[bin_id][0] = new_record;

    // ground truth
    match_info mi = {
                std::to_string(bin_id), // bin_id
                bins[bin_id][seq_idx].id(), // ref_id
                query_id, // query_id
                match_len, // match_len
                insert_start_pos, // start_pos_in_ref
                query_match_start_pos // start_pos_in_query
            };

    return mi;
}


void run_program(cmd_arguments const & args)
{
    uint32_t match_len = args.match_len / 3; // because proteins
    // random number generator setup
    std::random_device rd;
    std::mt19937 gen(rd());

    // IO
    seqan3::sequence_file_input<aa_input_trait> q_fin{args.query_file_path};
    seqan3::sequence_file_input<aa_input_trait> r_fin{args.ref_file_path};
    //args.ref_out_dir do stuff with that to write to output


    // create ref seq vector ordered by bin
    std::vector<std::vector<sequence_record_type>> bins;
    uint32_t idx{0};
    uint32_t bin_idx{0};
    uint32_t offset = (args.ref_num % args.bin_num) * ((args.ref_num / args.bin_num) +1);
    std::vector<sequence_record_type> current_records;
    for (auto & rec : r_fin)
    {
        sequence_record_type cur_rec{std::move(rec.sequence()), std::move(rec.id())};
        current_records.push_back(cur_rec);
        ++idx;
        if ((bin_idx < args.ref_num % args.bin_num && idx % ((args.ref_num / args.bin_num) +1) == 0) || 
        (!(bin_idx < args.ref_num % args.bin_num)) && (idx - offset) % (args.ref_num / args.bin_num) == 0)
        {
            bins.push_back(current_records);
            current_records.clear();
            ++bin_idx;
        }
    }

    // create random array from 0 to bin_num-1 to choose from
    std::vector<uint32_t> nums;
    for (uint32_t i=0; i<args.bin_num; ++i) {nums.push_back(i);}
    std::shuffle(nums.begin(), nums.end(), gen);

    // create ground truth vector
    std::vector<match_info> match_infos;

    // iterate through queries, generate and potentially insert matches
    idx = 0;
    for (auto const & [sequence, id, quality] : q_fin)
    {
        // generate match
        uint64_t const query_length = std::ranges::size(sequence);
        std::uniform_int_distribution<uint64_t> match_start_dis(0, query_length - match_len);
        uint64_t const match_start_pos = match_start_dis(gen);
        std::vector<seqan3::aa27> match = sequence |
                                          seqan3::views::slice(match_start_pos, match_start_pos + match_len) |
                                          seqan3::ranges::to<std::vector>();

        // introduce errors
        uint32_t max_errors = (int) (match_len * args.max_error_rate);
        std::uniform_int_distribution<uint32_t> match_error_pos_dis(0, match_len -1);
        std::uniform_int_distribution<uint8_t> aa27_rank_dis(0, 26);
        for (uint8_t error_count = 0; error_count < max_errors; ++error_count)
        {
            uint32_t const error_pos = match_error_pos_dis(gen);
            seqan3::aa27 const current_base = match[error_pos];
            seqan3::aa27 new_base = current_base;
            while (new_base == current_base)
                seqan3::assign_rank_to(aa27_rank_dis(gen), new_base);
            match[error_pos] = new_base;
        }

        // choose insert mode
        std::uniform_int_distribution<uint64_t> insert_mode_dis(0, 9);
        uint64_t const insert_mode = insert_mode_dis(gen);

        //insert
        uint64_t insert_idx{0};
        if (insert_mode < 4) {continue;}
        else if (insert_mode < 7)
        {
            //insert to 1 bin
            match_info mi = insert_match(bins, match, nums[idx], 0, match_len, id, match_start_pos);
            match_infos.push_back(mi);
            ++idx;
        }
        else if (insert_mode < 8)
        {
            //insert to 1 bin 2 seqs
            match_info mi = insert_match(bins, match, nums[idx], 0, match_len, id, match_start_pos);
            match_infos.push_back(mi);
            match_info mi2 = insert_match(bins, match, nums[idx], 1, match_len, id, match_start_pos);
            match_infos.push_back(mi2);
            ++idx;
        }
        else if (insert_mode < 9)
        {
            //insert twice to same seq
            std::uniform_int_distribution<uint64_t> insert_start_dis1(0, std::ranges::size(bins[idx][0].sequence())/2 - match_len);
            uint64_t insert_start_pos = insert_start_dis1(gen);
            uint64_t insert_start_pos2 = std::ranges::size(bins[idx][0].sequence())/2 + insert_start_pos;
            int bin_id = nums[idx];
            std::vector<seqan3::aa27> seq_pre = bins[bin_id][0].sequence() |
                                                seqan3::views::slice(0, insert_start_pos) | seqan3::ranges::to<std::vector>();
            std::vector<seqan3::aa27> seq_mid = bins[bin_id][0].sequence() |
                                                seqan3::views::slice(insert_start_pos + match_len, insert_start_pos2) |
                                                seqan3::ranges::to<std::vector>();
            std::vector<seqan3::aa27> seq_suf = bins[bin_id][0].sequence() |
                                                seqan3::views::slice(insert_start_pos2 + match_len, std::ranges::size(bins[bin_id][0].sequence())) | seqan3::ranges::to<std::vector>();
            std::vector<std::vector<seqan3::aa27>> v{ seq_pre, match, seq_mid, match, seq_suf };
            std::vector<seqan3::aa27> new_seq = std::ranges::join_view(v) |
                                                seqan3::ranges::to<std::vector>();

            // write new query to output file
            sequence_record_type new_record{std::move(new_seq), std::move(bins[bin_id][0].id())};
            bins[bin_id][0] = new_record;

            match_info mi = {
                std::to_string(bin_id), // bin_id
                bins[bin_id][0].id(), // ref_id
                id, // query_id
                match_len, // match_len
                insert_start_pos, // start_pos_in_ref
                match_start_pos // start_pos_in_query
            };
            match_info mi2 = {
                std::to_string(bin_id), // bin_id
                bins[bin_id][0].id(), // ref_id
                id, // query_id
                match_len, // match_len
                insert_start_pos2, // start_pos_in_ref
                match_start_pos // start_pos_in_query
            };
            match_infos.push_back(mi);
            match_infos.push_back(mi2);
            ++idx;
        }
        else if (insert_mode == 9)
        {
            //insert match to 2 bins 2 seqs
            match_info mi = insert_match(bins, match, nums[idx], 0, match_len, id, match_start_pos);
            match_infos.push_back(mi);
            ++idx;
            match_info mi2 = insert_match(bins, match, nums[idx], 0, match_len, id, match_start_pos);
            match_infos.push_back(mi2);
        }
    }

    // write bin-ref output files
    for (int i = 0; i < bins.size(); i++)
    {
        std::string filename = "bin_" + std::to_string(i) + ".fasta";
        std::filesystem::path out_file = args.ref_out_dir / filename;
        seqan3::sequence_file_output fout{out_file};
        for (auto & rec : bins[i])
        {
            fout.push_back(rec);
        }
    }
    
    // write ground truth file
    // sort by binId
    std::sort(match_infos.begin(), match_infos.end(), [](match_info a, match_info b)
                                                      {
                                                        return std::stoi(a.bin_id) < std::stoi(b.bin_id);
                                                      });

    std::ofstream ground_truth_file(args.ground_truth_file_path, std::ios_base::trunc);
    ground_truth_file << "bin_id, ref_id, query_id, match_len, start_pos_in_ref, start_pos_in_query \n";
    for (match_info info : match_infos)
    {
        ground_truth_file << info.bin_id << ", " << info.ref_id << ", " << info.query_id << ", " << info.match_len << ", " << info.start_pos_in_ref << ", " << info.start_pos_in_query << "\n";
    }
    
}

void initialize_argument_parser(seqan3::argument_parser & parser, cmd_arguments & args)
{
    // add meta data
    parser.info.app_name = "distribute_query_matches";
    parser.info.short_description = "distribute matches from queries into references";

    // add options/flags
    parser.add_option(args.query_file_path,
                      'q',
                      "in_queries",
                      "Please provide a file with query sequences.",
                      seqan3::option_spec::required,
                      seqan3::input_file_validator{{"fa", "fasta", "fna", "ffn", "ffa", "frn"}});

    parser.add_option(args.ref_file_path,
                      'r',
                      "in_refs",
                      "Please provide a file with reference sequences.",
                      seqan3::option_spec::required,
                      seqan3::input_file_validator{{"fa", "fasta", "fna", "ffn", "ffa", "frn"}});

    parser.add_option(args.ref_out_dir,
                      'o',
                      "out",
                      "Please provide an output directory in which the bins with reference sequences with inserted matches should be written.",
                      seqan3::option_spec::required, seqan3::output_directory_validator{});
                      
    parser.add_option(args.ground_truth_file_path,
                      'g',
                      "ground_truth",
                      "Please provide a file path for the ground truth file");

    parser.add_option(args.bin_num,
                      'b',
                      "bin_num",
                      "Please provide the number of bins.");

    parser.add_option(args.match_len,
                      'l',
                      "len_of_match",
                      "Please provide the length that the generated match in the query sequence should have.");

    parser.add_option(args.ref_num,
                      'n',
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
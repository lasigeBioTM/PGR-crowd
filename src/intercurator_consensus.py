import csv
from external_rater import n, k, N, join_amazon_external_rater_fleiss_kappa


def create_consensus_dataset(dataset_file, external_rater_file, destination_file):

    kappa, dict_counts = join_amazon_external_rater_fleiss_kappa(dataset_file, external_rater_file, n, N)

    dataset = open(dataset_file, encoding='utf-8')
    dataset_reader = csv.reader(dataset, delimiter=',')

    line_count = 0

    dict_content = {}

    for row in dataset_reader:

        if line_count == 0:
            header = row

        elif row:
            if row[21] == '' and row[27] not in dict_content:
                dict_content[row[27]] = row

        line_count += 1

    dataset.close()

    sentence_list = []
    blacklist = []
    for sentence, counts in dict_counts.items():
        for count in counts:
            if count in [6, 7, 8] and count == counts[0]:
                sentence_list.append(dict_content[sentence][:-1] + ['Yes, they share a direct/explicit relation in the sentence.'])
            elif count in [6, 7, 8] and count == counts[1]:
                sentence_list.append( dict_content[sentence][:-1] + ['No, they are separate entities with no correlation in the sentence.'])
            elif count in [6, 7, 8] and count == counts[2]:
                sentence_list.append(dict_content[sentence][:-1] + ['The entities seem to be illy marked, or something is wrong with the entities/sentence.'])
        if counts[2] == 4 or counts[2] == 5:
            blacklist.append(sentence)
        elif 6 not in counts and 7 not in counts and 8 not in counts and counts[2] != 4 and counts[2] != 5:
            sentence_list.append(dict_content[sentence][:-1] + ['No, they are separate entities with no correlation in the sentence.'])

    output_file = open(destination_file, 'w', encoding='utf-8')
    output_file_writer = csv.writer(output_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_ALL)

    output_file_writer.writerow(header)
    for row in sentence_list:
        output_file_writer.writerow(row)

    return blacklist


#create_consensus_dataset('data/batch_results_30.csv', 'data/external_rater_results.tsv', 'data/batch_results_30_consensus.csv')

import csv
import krippendorff
import numpy


### EXTERNAL VARIABLES

n = 2  # curators
N = 2389  # subjects
k = 3  # categories


def fleiss_kappa(dataset_file, curators, subjects):
    """

    :param dataset_file:
    :param curators:
    :param subjects:
    :return:
    """

    dataset = open(dataset_file, encoding = 'utf-8')
    dataset_reader = csv.reader(dataset, delimiter = ',')

    line_count = 0
    dict_counts = {}
    sentence_number = 1

    t_count = 0
    f_count = 0
    u_count = 0

    for row in dataset_reader:
        if line_count == 0:
            pass
        elif row and row[0] != '367O8HRHKG9KGF382XKL72J6UFOS44':
            if row[21] == '' and t_count + f_count + u_count < curators:  # not rejected and divisible by the number of curators:

                if row[28].startswith('Yes'):
                    t_count += 1
                elif row[28].startswith('No'):
                    f_count += 1
                elif row[28].startswith('The'):
                    u_count += 1

            elif row[21] == '' and t_count + f_count + u_count >= curators:

                dict_counts[sentence_number] = [t_count, f_count, u_count]
                sentence_number += 1
                t_count = 0
                f_count = 0
                u_count = 0

                if row[28].startswith('Yes'):
                    t_count += 1
                elif row[28].startswith('No'):
                    f_count += 1
                elif row[28].startswith('The'):
                    u_count += 1

        elif row and row[0] == '367O8HRHKG9KGF382XKL72J6UFOS44':
            if sentence_number not in dict_counts:
                dict_counts[sentence_number] = [7, 0, 0]
                sentence_number += 1

        line_count += 1

    sum_pi = 0
    t_sum = 0
    f_sum = 0
    u_sum = 0

    for item in dict_counts.items():
        pi = (1 / (curators * (curators - 1))) * (item[1][0]**2 + item[1][1]**2 + item[1][2]**2 - curators)
        sum_pi += pi
        t_sum += item[1][0]
        f_sum += item[1][1]
        u_sum += item[1][2]

    t = t_sum / (subjects * curators)
    f = f_sum / (subjects * curators)
    u = u_sum / (subjects * curators)

    p = (1 / subjects) * sum_pi
    pe = t**2 + f**2 + u**2

    kappa = (p - pe) / (1 - pe)

    return kappa

#print(fleiss_kappa('data/batch_results_30.csv', n, N))


def krippendorff_alpha(dataset_file, curators):
    """

    :param dataset_file:
    :param curators:
    :return:
    """

    dataset = open(dataset_file, encoding = 'utf-8')
    dataset_reader = csv.reader(dataset, delimiter = ',')

    line_count = 0
    reliability_data = []

    t_count = 0
    f_count = 0
    u_count = 0
    t = 1
    f = 2
    u = 3

    subset = []
    entered = False

    for row in dataset_reader:

        if line_count == 0:
            pass
        elif row and row[0] != '367O8HRHKG9KGF382XKL72J6UFOS44':
            if row[21] == '' and t_count + f_count + u_count < curators:  # not rejected and divisible by the number of curators:

                if row[28].startswith('Yes'):
                    subset.append(t)
                    t_count += 1
                elif row[28].startswith('No'):
                    subset.append(f)
                    f_count += 1
                elif row[28].startswith('The'):
                    subset.append(u)
                    u_count += 1

            elif row[21] == '' and t_count + f_count + u_count >= curators:

                reliability_data.append(subset)
                subset = []
                t_count = 0
                f_count = 0
                u_count = 0

                if row[28].startswith('Yes'):
                    subset.append(t)
                    t_count += 1
                elif row[28].startswith('No'):
                    subset.append(f)
                    f_count += 1
                elif row[28].startswith('The'):
                    subset.append(u)
                    u_count += 1

        elif row and row[0] == '367O8HRHKG9KGF382XKL72J6UFOS44' and entered == False:
            reliability_data.append([1, 1, 1, 1, 1, 1, 1])
            entered = True

        line_count += 1

    reliability_data = numpy.transpose(reliability_data)

    return krippendorff.alpha(reliability_data, level_of_measurement = 'nominal')

#print(krippendorff_alpha('data/batch_results_30.csv', n))

def expert_extra_kappa(external_rater_file, expert_rater_file, curators, subjects):

    dict_counts = {}
    t_count = 0
    f_count = 0
    u_count = 0

    external_rater = open(external_rater_file, encoding='utf-8')
    external_rater_reader = csv.reader(external_rater, delimiter='\t')

    line_count = 0

    for external_row in external_rater_reader:
        t_count = 0
        f_count = 0
        u_count = 0
        if line_count == 0:
            pass
        else:
            grading = external_row[11]
            gene_first = external_row[6]
            gene_last = external_row[7]
            phenotype_first = external_row[8]
            phenotype_last = external_row[9]

            if int(gene_first) < int(phenotype_first):
                new_sentence = external_row[1][:int(gene_first)] + '<b>' + external_row[1][
                                                                           int(gene_first):int(gene_last)] + '</b>' + \
                               external_row[1][int(gene_last):int(phenotype_first)] + '<b>' + external_row[1][
                                                                                              int(phenotype_first):int(
                                                                                                  phenotype_last)] + \
                               '</b>' + external_row[1][int(phenotype_last):]

            else:
                new_sentence = external_row[1][:int(phenotype_first)] + '<b>' + external_row[1][
                                                                                int(phenotype_first):int(
                                                                                    phenotype_last)] + '</b>' + \
                               external_row[1][int(phenotype_last):int(gene_first)] + '<b>' + external_row[1][
                                                                                              int(gene_first):int(
                                                                                                  gene_last)] + \
                               '</b>' + external_row[1][int(gene_last):]

            sentence = new_sentence.replace(',', '<span>&#44;</span>')

            if grading == 'C':
                t_count = 1
            elif grading == 'I':
                f_count = 1
            elif grading == 'U':
                u_count = 1
            dict_counts[sentence] = [t_count, f_count, u_count]

        line_count += 1


    expert_rater = open(expert_rater_file, encoding='utf-8')
    expert_rater_reader = csv.reader(expert_rater, delimiter='\t')

    line_count = 0

    for external_row in expert_rater_reader:
        if line_count == 0:
            pass
        else:
            grading = external_row[11].capitalize()
            gene_first = external_row[6]
            gene_last = external_row[7]
            phenotype_first = external_row[8]
            phenotype_last = external_row[9]

            if int(gene_first) < int(phenotype_first):
                new_sentence = external_row[1][:int(gene_first)] + '<b>' + external_row[1][
                                                                           int(gene_first):int(gene_last)] + '</b>' + \
                               external_row[1][int(gene_last):int(phenotype_first)] + '<b>' + external_row[1][
                                                                                              int(phenotype_first):int(
                                                                                                  phenotype_last)] + \
                               '</b>' + external_row[1][int(phenotype_last):]

            else:
                new_sentence = external_row[1][:int(phenotype_first)] + '<b>' + external_row[1][
                                                                                int(phenotype_first):int(
                                                                                    phenotype_last)] + '</b>' + \
                               external_row[1][int(phenotype_last):int(gene_first)] + '<b>' + external_row[1][
                                                                                              int(gene_first):int(
                                                                                                  gene_last)] + \
                               '</b>' + external_row[1][int(gene_last):]

            sentence = new_sentence.replace(',', '<span>&#44;</span>')

            if grading == 'C':
                dict_counts[sentence] = [dict_counts[sentence][0] + 1, dict_counts[sentence][1], dict_counts[sentence][2]]
            elif grading == 'I':
                dict_counts[sentence] = [dict_counts[sentence][0], dict_counts[sentence][1] + 1, dict_counts[sentence][2]]
            elif grading.startswith('U'):
                dict_counts[sentence] = [dict_counts[sentence][0], dict_counts[sentence][1], dict_counts[sentence][2] + 1]

        line_count += 1

    sum_pi = 0
    t_sum = 0
    f_sum = 0
    u_sum = 0

    for item in dict_counts.items():
        # if item[1][0] + item[1][1] + item[1][2] != 2:
        #     print(item)
        # elif item[1][1] == 1 and item[1][2] == 1:
        #     print(item)
        pi = (1 / (curators * (curators - 1))) * (item[1][0] ** 2 + item[1][1] ** 2 + item[1][2] ** 2 - curators)
        sum_pi += pi
        t_sum += item[1][0]
        f_sum += item[1][1]
        u_sum += item[1][2]

    t = t_sum / (subjects * curators)
    f = f_sum / (subjects * curators)
    u = u_sum / (subjects * curators)

    p = (1 / subjects) * sum_pi
    pe = t ** 2 + f ** 2 + u ** 2

    kappa = (p - pe) / (1 - pe)

    return kappa

print(expert_extra_kappa('data/external_rater_results.tsv', 'data/expert.tsv', n, N))


def expert_extra_krip(external_rater_file, expert_rater_file):

    data = {}

    t = 1
    f = 2
    u = 3

    expert_rater = open(expert_rater_file, encoding='utf-8')
    expert_rater_reader = csv.reader(expert_rater, delimiter='\t')

    line_count = 0

    for external_row in expert_rater_reader:
        if line_count == 0:
            pass
        else:
            grading = external_row[11].capitalize()
            gene_first = external_row[6]
            gene_last = external_row[7]
            phenotype_first = external_row[8]
            phenotype_last = external_row[9]

            if int(gene_first) < int(phenotype_first):
                new_sentence = external_row[1][:int(gene_first)] + '<b>' + external_row[1][
                                                                           int(gene_first):int(gene_last)] + '</b>' + \
                               external_row[1][int(gene_last):int(phenotype_first)] + '<b>' + external_row[1][
                                                                                              int(phenotype_first):int(
                                                                                                  phenotype_last)] + \
                               '</b>' + external_row[1][int(phenotype_last):]

            else:
                new_sentence = external_row[1][:int(phenotype_first)] + '<b>' + external_row[1][
                                                                                int(phenotype_first):int(
                                                                                    phenotype_last)] + '</b>' + \
                               external_row[1][int(phenotype_last):int(gene_first)] + '<b>' + external_row[1][
                                                                                              int(gene_first):int(
                                                                                                  gene_last)] + \
                               '</b>' + external_row[1][int(gene_last):]

            sentence = new_sentence.replace(',', '<span>&#44;</span>')
            if grading == 'C':
                data[sentence] = [t]
            elif grading == 'I':
                data[sentence] = [f]
            elif grading.startswith('U'):
                data[sentence] = [u]

        line_count += 1

    external_rater = open(external_rater_file, encoding='utf-8')
    external_rater_reader = csv.reader(external_rater, delimiter='\t')

    line_count = 0

    for external_row in external_rater_reader:
        if line_count == 0:
            pass
        else:
            grading = external_row[11]
            gene_first = external_row[6]
            gene_last = external_row[7]
            phenotype_first = external_row[8]
            phenotype_last = external_row[9]

            if int(gene_first) < int(phenotype_first):
                new_sentence = external_row[1][:int(gene_first)] + '<b>' + external_row[1][
                                                                           int(gene_first):int(gene_last)] + '</b>' + \
                               external_row[1][int(gene_last):int(phenotype_first)] + '<b>' + external_row[1][
                                                                                              int(phenotype_first):int(
                                                                                                  phenotype_last)] + \
                               '</b>' + external_row[1][int(phenotype_last):]

            else:
                new_sentence = external_row[1][:int(phenotype_first)] + '<b>' + external_row[1][
                                                                                int(phenotype_first):int(
                                                                                    phenotype_last)] + '</b>' + \
                               external_row[1][int(phenotype_last):int(gene_first)] + '<b>' + external_row[1][
                                                                                              int(gene_first):int(
                                                                                                  gene_last)] + \
                               '</b>' + external_row[1][int(gene_last):]

            sentence = new_sentence.replace(',', '<span>&#44;</span>')
            if grading == 'C':
                data[sentence].append(t)
            elif grading == 'I':
                data[sentence].append(f)
            elif grading == 'U':
                data[sentence].append(u)

        line_count += 1

    reliability_data = []
    for item in data.items():
        reliability_data.append(item[1])

    reliability_data = numpy.transpose(reliability_data)

    external_rater.close()

    return krippendorff.alpha(reliability_data, level_of_measurement='nominal')

print(expert_extra_krip('data/external_rater_results.tsv', 'data/expert.tsv'))
import csv
import krippendorff
import numpy


### EXTERNAL VARIABLES

n = 7  # curators
N = 2389  # subjects
k = 3  # categories


def join_amazon_external_rater_fleiss_kappa(dataset_file, external_rater_file, curators, subjects):
    """

    :param dataset_file:
    :param external_rater_file:
    :param curators:
    :param subjects
    :return:
    """

    dataset = open(dataset_file, encoding='utf-8')
    dataset_reader = csv.reader(dataset, delimiter=',')

    line_count = 0
    dict_counts = {}

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

                if t_count + f_count + u_count == curators:
                    sentence = row[27]
                    dict_counts[sentence] = [t_count, f_count, u_count]

            elif row[21] == '' and t_count + f_count + u_count >= curators:

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
            if row[27] not in dict_counts:
                dict_counts[row[27]] = [7, 0, 0]

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
                new_sentence = external_row[1][:int(gene_first)] + '<b>' + external_row[1][int(gene_first):int(gene_last)] + '</b>' + \
                               external_row[1][int(gene_last):int(phenotype_first)] + '<b>' + external_row[1][int(phenotype_first):int(phenotype_last)] + \
                               '</b>' + external_row[1][int(phenotype_last):]

            else:
                new_sentence = external_row[1][:int(phenotype_first)] + '<b>' + external_row[1][int(phenotype_first):int(phenotype_last)] + '</b>' + \
                               external_row[1][int(phenotype_last):int(gene_first)] + '<b>' + external_row[1][int(gene_first):int(gene_last)] + \
                               '</b>' + external_row[1][int(gene_last):]

            sentence = new_sentence.replace(',', '<span>&#44;</span>')
            if grading == 'C':
                t_count = dict_counts[sentence][0] + 1
                dict_counts[sentence] = [t_count, dict_counts[sentence][1], dict_counts[sentence][2]]
            elif grading == 'I':
                f_count = dict_counts[sentence][1] + 1
                dict_counts[sentence] = [dict_counts[sentence][0], f_count, dict_counts[sentence][2]]
            elif grading == 'U':
                u_count = dict_counts[sentence][2] + 1
                dict_counts[sentence] = [dict_counts[sentence][0], dict_counts[sentence][1], u_count]

        line_count += 1

    sum_pi = 0
    t_sum = 0
    f_sum = 0
    u_sum = 0
    curators = 8

    for item in dict_counts.items():
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

    dataset.close()
    external_rater.close()

    return kappa, dict_counts

#kappa, dict_counts = join_amazon_external_rater_fleiss_kappa('data/batch_results_30.csv', 'data/external_rater_results.tsv', n, N)
#print(kappa)


def join_amazon_external_rater_krippendorff_alpha(dataset_file, external_rater_file, curators):
    """

    :param dataset_file:
    :param: external_rater_file:
    :param curators:
    :return:
    """

    dataset = open(dataset_file, encoding='utf-8')
    dataset_reader = csv.reader(dataset, delimiter=',')

    line_count = 0
    data = {}

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

                if t_count + f_count + u_count == curators:
                    sentence = row[27]
                    data[sentence] = subset

            elif row[21] == '' and t_count + f_count + u_count >= curators:

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

        elif row and row[0] == '367O8HRHKG9KGF382XKL72J6UFOS44' and entered==False:
            data[row[27]] = [1, 1, 1, 1, 1, 1, 1]
            entered = True

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
                new_sentence = external_row[1][:int(gene_first)] + '<b>' + external_row[1][ int(gene_first):int(gene_last)] + '</b>' + \
                               external_row[1][int(gene_last):int(phenotype_first)] + '<b>' + external_row[1][int(phenotype_first):int(phenotype_last)] + \
                               '</b>' + external_row[1][int(phenotype_last):]

            else:
                new_sentence = external_row[1][:int(phenotype_first)] + '<b>' + external_row[1][int(phenotype_first):int(phenotype_last)] + '</b>' + \
                               external_row[1][int(phenotype_last):int(gene_first)] + '<b>' + external_row[1][int(gene_first):int(gene_last)] + \
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

    dataset.close()
    external_rater.close()

    return krippendorff.alpha(reliability_data, level_of_measurement='nominal')

#print(join_amazon_external_rater_krippendorff_alpha('data/batch_results_30.csv', 'data/external_rater_results.tsv', n))


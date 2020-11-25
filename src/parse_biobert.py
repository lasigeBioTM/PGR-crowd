import os
import csv
import random
import numpy
import xml.etree.ElementTree as ET


def parse_tsv(base_dir, destination_file, nomenclature='relation', test=None):
    """

    :param base_dir:
    :param destination_file:
    :param nomenclature:
    :param test
    :return:
    """

    count_wrong = 0
    line_saves = []
    sentence_number = 0

    for f in os.listdir(base_dir):
        reader = open(base_dir + '/' + f, 'r', encoding='utf-8')
        content = reader.read().replace('</sup>', 'AAAA').replace('<sup>', 'AAA').replace('</b>', 'AAAA') \
            .replace('<b>', 'AAA').replace('</i>', 'AAAA').replace('<i>', 'AAA').replace('</sub>', 'AAAA') \
            .replace('<sub>', 'AAA')

        root = ET.fromstring(content)

        # tree = ET.parse(base_dir + '/' + f)
        # root = tree.getroot()

        for sentence in root:
            sentence_text = sentence.get('text')
            all_pairs = sentence.findall('pair')

            for pair in all_pairs:

                id_1 = pair.get('e1')
                id_2 = pair.get('e2')

                char_offset_1 = None
                char_offset_2 = None
                e_type_1 = None
                e_type_2 = None

                for e in sentence.findall('entity'):

                    if e.get('id') == id_1:
                        if e.get('type') == 'GO':
                            e_type_1 = 'GENE'
                        else:
                            e_type_1 = 'DISEASE'

                        char_offset_1 = e.get('charOffset')

                    elif e.get('id') == id_2:
                        if e.get('type') == 'GO':
                            e_type_2 = 'GENE'
                        else:
                            e_type_2 = 'DISEASE'

                        char_offset_2 = e.get('charOffset')

                sentence_final = sentence_text[:int(char_offset_1.split('-')[0])] + '@' + e_type_1 + '$'\
                               + sentence_text[int(char_offset_1.split('-')[1]):int(char_offset_2.split('-')[0])]\
                               + '@' + e_type_2 + '$' + sentence_text[int(char_offset_2.split('-')[1]):]

                relation = pair.get(nomenclature)
                if relation == 'true':
                    label = 1
                else:
                    label = 0

                if test:
                    line_saves.append([sentence_number, sentence_final, label])
                else:
                    line_saves.append([sentence_final, label])
                sentence_number += 1

    output_file = open(destination_file + '.tsv', 'w', encoding='utf-8')
    output_file_writer = csv.writer(output_file, delimiter='\t')

    if test:
        output_file_writer.writerow(['index', 'sentence', 'label'])

        for row in line_saves:
            output_file_writer.writerow(row)

    else:
        dev_output_file = open(destination_file + '_dev.tsv', 'w', encoding='utf-8')
        dev_output_file_writer = csv.writer(dev_output_file, delimiter='\t')

        dev = random.sample(line_saves, int((len(line_saves) * 20) / 100))
        for row in dev:
            dev_output_file_writer.writerow(row)
        for row in line_saves:
            if row not in dev:
                output_file_writer.writerow(row)

    return


#parse_tsv('corpora/amazon_train', 'corpora/amazon_train', nomenclature='relation')
#parse_tsv('corpora/consensus_test', 'corpora/consensus_train', nomenclature='relation')
#parse_tsv('corpora/consensus_test', 'corpora/consensus_test', nomenclature='relation', test=True)
#parse_tsv('corpora/original_test', 'corpora/original_test', nomenclature='relation', test=True)
#parse_tsv('corpora/original_train', 'corpora/original_train', nomenclature='relation')
#parse_tsv('corpora/expert_test', 'corpora/expert_test', nomenclature='relation', test=True)
#parse_tsv('corpora/expert_test', 'corpora/expert_train', nomenclature='relation')

# import pandas  as pd
# import matplotlib.pyplot as plt
# import numpy as np
# from sklearn.linear_model import LogisticRegression
# from sklearn.preprocessing import StandardScaler, LabelEncoder
# from sklearn.metrics import confusion_matrix, classification_report
# from sklearn.model_selection import train_test_split
# from imblearn.over_sampling import SMOTE
#
# def smote(dataset_file):
#
#     data = pd.read_csv(dataset_file, sep='\t', names=['sentence', 'label'])
#     #print(data.info())
#
#     data['normalize_sentence'] =  LabelEncoder().fit_transform(np.array(data['sentence']).reshape(-1, 1).ravel())
#     print(data['label'].value_counts())
#
#     X_train = np.array(data['normalize_sentence'])
#     y_train = np.array(data['label'])
#
#     # describes info about train and test set
#     print("Number transactions X_train dataset: ", X_train.shape)
#     print("Number transactions y_train dataset: ", y_train.shape)
#
#     print("Before OverSampling, counts of label '1': {}".format(sum(y_train == 1)))
#     print("Before OverSampling, counts of label '0': {} \n".format(sum(y_train == 0)))
#
#     sm = SMOTE(random_state=2)
#     X_train_res, y_train_res = sm.fit_sample(X_train.reshape(-1, 1), y_train.ravel())
#
#     print('After OverSampling, the shape of train_X: {}'.format(X_train_res.shape))
#     print('After OverSampling, the shape of train_y: {} \n'.format(y_train_res.shape))
#
#     print("After OverSampling, counts of label '1': {}".format(sum(y_train_res == 1)))
#     print("After OverSampling, counts of label '0': {}".format(sum(y_train_res == 0)))
#
#     return
#
# smote('corpora/amazon_train.tsv')
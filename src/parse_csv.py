import csv
import os
import time
import collections
import xml.etree.ElementTree as ET
from xml.dom import minidom
from multiprocessing import Process


from intercurator_consensus import create_consensus_dataset, n, k, N




def get_70_dataset(original_dataset_file, original_dataset_30_file):
    """

    :param original_dataset_file:
    :param original_dataset_30_file:
    """

    dataset = open(original_dataset_file, encoding='utf-8')
    dataset_reader = csv.reader(dataset, delimiter='\t')

    all_sentences = []

    for row in dataset_reader:
        all_sentences.append(row)

    dataset.close()

    dataset_30 = open(original_dataset_30_file, encoding='utf-8')
    dataset_30_reader = csv.reader(dataset_30, delimiter='\t')

    line_count = 0
    sentences_30 = []

    for row in dataset_30_reader:
        if line_count == 0:
            pass
        else:
            sentences_30.append(row)

        line_count += 1

    dataset_30.close()

    dataset_70 = open('data/original_dataset_70.tsv', 'w', encoding='utf-8')

    for sentence in all_sentences:
        if sentence not in sentences_30:
            dataset_70.write('\t'.join(sentence) + '\n')

    dataset_70.close()

    return


get_70_dataset('data/original_dataset.tsv', 'data/original_dataset_30.tsv')


# GENE - GO DICTIONARY: GENE ID - (GO ID, EVIDENCE CODE, GO NAME, CATEGORY), (...) ####

def dict_g2go(file_g2go):
    """Creates a dictionary of type {gene1:[(GO_ID, Evidence, GO_name, category),
    (GO_ID, Evidence, GO_name, category), ...], }

    :param file_g2go: file with relations gene to GO
    :return: dict of type {gene1:[(GO_ID, Evidence, GO_name, category),
             (GO_ID, Evidence, GO_name, category), ...], }
    """

    os.system('gunzip -k ' + file_g2go + '.gz')

    gene2go = open(file_g2go, 'r', encoding = 'utf-8')

    gene2go.readline()  # skip header

    relations_g2go = gene2go.readlines()
    gene2go.close()

    relations_g2go.pop()

    dict_gene_go = {}

    for line in relations_g2go:

        line = line.split('\t')

        gene_id = line[1]
        go = line[2]
        evidence = line[3]
        name = line[5]
        category = line[7][:-1]

        if gene_id not in dict_gene_go:
            dict_gene_go[gene_id] = []
            dict_gene_go[gene_id].append((go, evidence, name, category))

        else:
            dict_gene_go[gene_id].append((go, evidence, name, category))

    os.system('rm ' + file_g2go)

    return dict_gene_go


# GENE ANNOTATIONS 2 GO ANNOTATIONS ####

def go_annotations(file_g2go):
    """Creates a dictionary of type {gene_id:go_id, gene_id2:go_id, } and
       a dictionary of type {gene_name:go_name, gene_name:go_name, }

    :param file_g2go: file with relations gene to GO
    :return: a dictionary of type {gene_id:go_id, gene_id2:go_id, } and
             a dictionary of type {gene_id:go_name, gene_id2:go_name, }
    """

    dict_gene_id_go = dict_g2go(file_g2go)
    dict_gene_go_id = {}
    dict_gene_go_name = {}

    for gene, list_evidence in dict_gene_id_go.items():

        criteria_order = ['EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP', 'HTP',
                          'HDA', 'HMP', 'HGI', 'HEP', 'ISS', 'ISO', 'ISA',
                          'ISM', 'IGC', 'IBA', 'IBD', 'IKR', 'IRD', 'RCA',
                          'TAS', 'NAS', 'IC', 'ND', 'IEA']  # order criteria

        bins = collections.defaultdict(list)
        for pair in list_evidence:
            bins[pair[1]].append(pair)

        value_tup = [pair for i in criteria_order for pair in bins[i]]

        for v_tup in value_tup:

            if v_tup[3] == 'Process':  # Biological Process

                dict_gene_go_id[gene] = v_tup[0].replace(':', '_')
                dict_gene_go_name[gene] = v_tup[2]
                break

    return dict_gene_go_id, dict_gene_go_name


def get_pubmed_id_sentences(dataset_file):
    """

    :param dataset_file:
    :return:
    """

    dataset = open(dataset_file, encoding='utf-8')
    dataset_reader = csv.reader(dataset, delimiter='\t')

    line_count = 0
    dict_sentence_id = {}
    dict_id_sentence = {}

    for row in dataset_reader:
        if line_count == 0:
            pass
        else:
            dict_sentence_id[row[1]] = row[0]

            if row[0] in dict_id_sentence:
                if row[1] not in dict_id_sentence[row[0]]:
                    dict_id_sentence[row[0]].append(row[1])
            else:
                dict_id_sentence[row[0]] = [row[1]]

        line_count += 1

    dataset.close()

    return dict_id_sentence


def get_entities_sentence(dataset_file, sentence):
    """

    :param dataset_file:
    :param sentence:
    :return: sentence : [[entity1, id, char1, char2], [entity2, id, char3, char4],
                        relation], [[entity1, id, char1, char2], [entity2, id, char3, char4],
                        relation], etc.]
    """

    dataset = open(dataset_file, encoding='utf-8')
    dataset_reader = csv.reader(dataset, delimiter='\t')

    line_count = 0
    dict_entities_sentence = {}

    for row in dataset_reader:
        if line_count == 0:
            pass
        else:
            if row[1] in dict_entities_sentence:
                if row[6] > row[8]:
                    dict_entities_sentence[row[1]].append([[row[3], row[5], row[8], row[9]], [row[2], row[4], row[6], row[7]], row[10]])
                else:
                    dict_entities_sentence[row[1]].append([[row[2], row[4], row[6], row[7]], [row[3], row[5], row[8], row[9]], row[10]])
            else:

                if row[6] > row[8]:
                    dict_entities_sentence[row[1]] = [[[row[3], row[5], row[8], row[9]], [row[2], row[4], row[6], row[7]], row[10]]]
                else:
                    dict_entities_sentence[row[1]] = [[[row[2], row[4], row[6], row[7]], [row[3], row[5], row[8], row[9]], row[10]]]

        line_count += 1

    dataset.close()

    return dict_entities_sentence[sentence]

def get_pubmed_id_sentences_expert(dataset_file):
    """

    :param dataset_file:
    :return:
    """

    dataset = open(dataset_file, encoding='utf-8')
    dataset_reader = csv.reader(dataset, delimiter='\t')

    line_count = 0
    dict_sentence_id = {}
    dict_id_sentence = {}

    for row in dataset_reader:
        if line_count == 0:
            pass
        else:
            if row[11] == 'I' or row[11] == 'C' or row[11] == 'i' or row[11] == 'c':
                dict_sentence_id[row[1]] = row[0]

                if row[0] in dict_id_sentence:
                    if row[1] not in dict_id_sentence[row[0]]:
                        dict_id_sentence[row[0]].append(row[1])
                else:
                    dict_id_sentence[row[0]] = [row[1]]

        line_count += 1

    dataset.close()

    return dict_id_sentence

def get_entities_sentence_expert(dataset_file, sentence):
    """

    :param dataset_file:
    :param sentence:
    :return: sentence : [[entity1, id, char1, char2], [entity2, id, char3, char4],
                        relation], [[entity1, id, char1, char2], [entity2, id, char3, char4],
                        relation], etc.]
    """

    dataset = open(dataset_file, encoding='utf-8')
    dataset_reader = csv.reader(dataset, delimiter='\t')

    line_count = 0
    dict_entities_sentence = {}

    for row in dataset_reader:
        if line_count == 0:
            pass
        else:
            if row[1] in dict_entities_sentence:
                relation = ''
                if (row[11] == 'C' or row[11] == 'c') and row[10] == '1':
                    relation = 'True'
                elif (row[11] == 'I' or row[11] == 'i') and row[10] == '0':
                    relation = 'True'
                elif (row[11] == 'C' or row[11] == 'c') and row[10] == '0':
                    relation = 'False'
                elif (row[11] == 'I' or row[11] == 'i') and row[10] == '1':
                    relation = 'False'
                elif row[11].startswith('U') or row[11].startswith('u'):
                    continue
                if row[6] > row[8]:
                    dict_entities_sentence[row[1]].append([[row[3], row[5], row[8], row[9]], [row[2], row[4], row[6], row[7]], relation])
                else:
                    dict_entities_sentence[row[1]].append([[row[2], row[4], row[6], row[7]], [row[3], row[5], row[8], row[9]], relation])
            else:
                relation = ''
                if (row[11] == 'C' or row[11] == 'c') and row[10] == '1':
                    relation = 'True'
                elif (row[11] == 'I' or row[11] == 'i') and row[10] == '0':
                    relation = 'True'
                elif (row[11] == 'C' or row[11] == 'c') and row[10] == '0':
                    relation = 'False'
                elif (row[11] == 'I' or row[11] == 'i') and row[10] == '1':
                    relation = 'False'
                elif row[11].startswith('U') or row[11].startswith('u'):
                    continue
                if row[6] > row[8]:
                    dict_entities_sentence[row[1]] = [[[row[3], row[5], row[8], row[9]], [row[2], row[4], row[6], row[7]], relation]]
                else:
                    dict_entities_sentence[row[1]] = [[[row[2], row[4], row[6], row[7]], [row[3], row[5], row[8], row[9]], relation]]

        line_count += 1

    dataset.close()

    return dict_entities_sentence[sentence]
# CREATE XML FILES ####

def prettify(elem):
    """Return a pretty-printed XML string for the Element

    :param elem:
    :return:
    """

    rough_string = ET.tostring(elem, 'utf-8')
    reparsed = minidom.parseString(rough_string)

    return reparsed.toprettyxml(indent='  ')


def xml_file(dataset_file, amazon_file, validation_file, file_g2go, destination_path, test=False):
    """Process to create each file

    :param dataset_file:
    :param amazon_file:
    :param validation_file:
    :param file_g2go:
    :param destination_path:
    :param test:
    """

    if test is True:
        blacklist = create_consensus_dataset('data/batch_results_30.csv', 'data/external_rater_results.tsv', 'data/batch_results_30_consensus.csv')
    else:
        blacklist = []

    # stats
    count = 0
    count_true = 0
    count_false = 0
    count_excluded = 0

    original_count = 0
    original_count_true = 0
    original_count_false = 0

    validation = open(validation_file, 'r', encoding='utf-8')
    validation.readline()
    validation_sentences = validation.readlines()

    dict_gene_go_id, dict_gene_go_name = go_annotations(file_g2go)
    amazon = open(amazon_file, encoding='utf-8')
    amazon_file_reader = csv.reader(amazon, delimiter=',', quotechar='"')

    line_count = 0
    dict_replies = {}  # sentence : reply

    for row in amazon_file_reader:
        if line_count == 0:
            pass
        elif row:
            if row[21] == '' and test is False and row[-2] + '\n' not in validation_sentences:
                if row[-1] == 'Yes, they share a direct/explicit relation in the sentence.':
                    dict_replies[row[-2]] = 'true'
                    count_true += 1
                elif row[-1] == 'The entities seem to be illy marked, or something is wrong with the entities/sentence.':
                    dict_replies[row[-2]] = 'error'
                    count_excluded += 1
                elif row[-1] == 'No, they are separate entities with no correlation in the sentence.':
                    dict_replies[row[-2]] = 'false'
                    count_false += 1
            elif row[21] == '' and test is True:
                if row[-1] == 'Yes, they share a direct/explicit relation in the sentence.':
                    dict_replies[row[-2]] = 'true'
                    count_true += 1
                elif row[-1] == 'The entities seem to be illy marked, or something is wrong with the entities/sentence.':
                    dict_replies[row[-2]] = 'error'
                    count_excluded += 1
                elif row[-1] == 'No, they are separate entities with no correlation in the sentence.':
                    dict_replies[row[-2]] = 'false'
                    count_false += 1

        line_count += 1

    amazon.close()

    dict_id_sentence = get_pubmed_id_sentences(dataset_file)  # id_pubmed : [sentence1, sentence2, etc.]

    for key, items in dict_id_sentence.items():
        root = ET.Element('document', id=key)
        sentence_number = 0

        for item in items:
            entity_number = 0
            entities_sentence = get_entities_sentence(dataset_file, item)

            doc = ET.SubElement(root, 'sentence', id=key + '.s' + str(sentence_number), text=item)
            save_entities = []
            save_pairs = []

            for pair in entities_sentence:  # sentence : [[entity1, id, char1, char2], [entity2, id, char3, char4],
                                            # relation], [[entity1, id, char1, char2], [entity2, id, char3, char4],
                                            # relation], etc.]

                pair_list = []

                entity_1 = pair[0]
                entity_2 = pair[1]

                if entity_1 not in save_entities:
                    save_entities.append(entity_1)

                if entity_2 not in save_entities:
                    save_entities.append(entity_2)

                pair_list.append(entity_1)
                pair_list.append(entity_2)

                pair_list_sorted = sorted(pair_list, key=lambda x: int(x[2]))
                pair_list_sorted.append(pair[2])

                save_pairs.append(pair_list_sorted)

            save_entities_sorted = sorted(save_entities, key = lambda x: int(x[2]))
            associated_entity_number = {}

            for entity in save_entities_sorted:

                if entity[1].startswith('HP'):

                    ET.SubElement(doc, 'entity', id=key + '.s' + str(sentence_number) + '.e' + str(entity_number),
                                  charOffset=entity[2] + '-' + entity[3], type='HP', text=entity[0], ontology_id=entity[1])

                else:

                    ET.SubElement(doc, 'entity', id=key + '.s' + str(sentence_number) + '.e' + str(entity_number),
                                  charOffset=entity[2] + '-' + entity[3], type='GENE', text=entity[0], ontology_id=entity[1])

                associated_entity_number[repr(entity)] = key + '.s' + str(sentence_number) + '.e' + str(entity_number)
                entity_number += 1

            # save_pairs = [[['PDYN', '5173', '64', '68'], ['epilepsy', 'HP_0001250', '111', '119'], 'False'], [['PDYN', '5173', '32', '68'], ['epilepsy', 'HP_0001250', '111', '119'], 'False']]

            pairs_sorted = sorted(save_pairs, key = lambda x: int(x[0][2]))
            pair_number = 0

            for pair in pairs_sorted:

                original_count += 1
                if pair[2] == 'False':
                    original_count_false += 1
                elif pair[2] == 'True':
                    original_count_true += 1

                new_sentence = item[:int(pair[0][2])] + '<b>' + item[int(pair[0][2]):int(pair[0][3])] + '</b>' + \
                               item[int(pair[0][3]):int(pair[1][2])] + '<b>' + item[int(pair[1][2]):int(pair[1][3])] + \
                               '</b>' + item[int(pair[1][3]):]

                if new_sentence.replace(',', '<span>&#44;</span>') in blacklist:
                    pass
                else:
                    if dict_replies[new_sentence.replace(',', '<span>&#44;</span>')] == 'error':
                         pass
                    else:
                        ET.SubElement(doc, 'pair', id=key + '.s' + str(sentence_number) + '.p' + str(pair_number),
                                e1=associated_entity_number[repr(pair[0])], e2=associated_entity_number[repr(pair[1])],
                                relation=dict_replies[new_sentence.replace(',', '<span>&#44;</span>')])
                        count += 1
                        pair_number += 1

            sentence_number += 1

        output_file = open(destination_path + key + '.xml', 'w', encoding='utf-8')
        # output_file.write('<?xml version="1.0" encoding="UTF-8"?>')

        root = xml_file_go(dict_gene_go_name, dict_gene_go_id, root)
        output_file.write(prettify(root))

        output_file.close()

    # stats
    print('total', count)
    print('count true', count_true)
    print('count false', count_false)
    print('count_excluded', count_excluded)

    print('------------------------------------------------------------')

    print('original count', original_count)
    print('original_count_true', original_count_true)
    print('original_count_false', original_count_false)

    return

def xml_file_go(dict_gene_go_name, dict_gene_go_id, xml_root):
    """

    :param dict_gene_go_name:
    :param dict_gene_go_id:
    :param xml_root:
    :return:
    """

    root = xml_root

    for sentence in root:
        save_genes = {}  # gene : [gene_id, char_1, char_2]
        save_added_length = 0

        for sub_element in sentence:
            if sub_element.get('type') == 'GENE':
                save_genes[sub_element.get('text')] = [sub_element.get('ontology_id'), sub_element.get('charOffset').split('-')[0], sub_element.get('charOffset').split('-')[1]]
                save_add = len(sub_element.get('text'))

                if sub_element.get('ontology_id') not in dict_gene_go_name:
                    sub_element.attrib['ontology_id'] = sub_element.attrib['ontology_id'].replace(sub_element.get('ontology_id'), 'GO_0008150')
                    sub_element.attrib['text'] = sub_element.attrib['text'].replace(sub_element.get('text'), 'biological_process')
                else:
                    save = sub_element.get('ontology_id')
                    sub_element.attrib['ontology_id'] = sub_element.attrib['ontology_id'].replace(sub_element.get('ontology_id'), dict_gene_go_id[sub_element.get('ontology_id')])
                    sub_element.attrib['text'] = sub_element.attrib['text'].replace(sub_element.get('text'), dict_gene_go_name[save])

                sub_element.attrib['type'] = sub_element.attrib['type'].replace(sub_element.get('type'), 'GO')
                sub_element.attrib['charOffset'] = sub_element.attrib['charOffset'].replace(sub_element.get('charOffset'), str(int(sub_element.get('charOffset').split('-')[0]) + save_added_length) + '-' + str(int(sub_element.get('charOffset').split('-')[0]) + save_added_length + len(sub_element.get('text'))))
                save_added_length -= save_add
                save_added_length += len(sub_element.get('text'))


            elif 'e' in sub_element.get('id'):
                sub_element.attrib['charOffset'] = sub_element.attrib['charOffset'].replace(sub_element.get('charOffset'), str(int(sub_element.get('charOffset').split('-')[0]) + save_added_length) + '-' + str(int(sub_element.get('charOffset').split('-')[1]) + save_added_length))

        save_added_length_sentence = 0
        for gene, gene_elements in save_genes.items():
            if gene_elements[0] not in dict_gene_go_name:
                sentence.attrib['text'] = sentence.attrib['text'].replace(sentence.get('text'), sentence.get('text')[:int(gene_elements[1]) + save_added_length_sentence]) + 'biological_process' + sentence.get('text')[int(gene_elements[2]) + save_added_length_sentence:]
                save_added_length_sentence -= len(gene)
                save_added_length_sentence += len('biological_process')
            else:
                sentence.attrib['text'] = sentence.attrib['text'].replace(sentence.get('text'), sentence.get('text')[:int(gene_elements[1]) + save_added_length_sentence]) + dict_gene_go_name[gene_elements[0]] + sentence.get('text')[int(gene_elements[2]) + save_added_length_sentence:]
                save_added_length_sentence -= len(gene)
                save_added_length_sentence += len(dict_gene_go_name[gene_elements[0]])

    return root

def get_original_test_set_xml(original_dataset_file, destination_file, destination_path, file_g2go):
    """

    :param original_dataset_file:
    :param destination_file:
    :param destination_path:
    :param file_g2go:
    :return:
    """
    dict_gene_go_id, dict_gene_go_name = go_annotations(file_g2go)

    # for test
    # dataset = open(original_dataset_file, encoding='utf-8')
    # dataset_reader = csv.reader(dataset, delimiter='\t')
    #
    # line_count = 0
    # line_saves = []
    #
    # for row in dataset_reader:
    #     if line_count == 0:
    #         line_saves.append(row[:-1])
    #     else:
    #         if row[11] == 'C':
    #             line_saves.append(row[:-1])
    #         elif row[11] == 'I' and row[10] == 'TRUE':
    #             line_saves.append(row[:-2] + ['FALSE'])
    #         elif row[11] == 'I' and row[10] == 'FALSE':
    #             line_saves.append(row[:-2] + ['TRUE'])
    #
    #     line_count += 1
    #
    # output_file = open(destination_file, 'w', encoding='utf-8')
    # output_file_writer = csv.writer(output_file, delimiter='\t')
    #
    # for row in line_saves:
    #     output_file_writer.writerow(row)
    # end for test

    #dict_id_sentence = get_pubmed_id_sentences(destination_file)  # id_pubmed : [sentence1, sentence2, etc.]
    dict_id_sentence = get_pubmed_id_sentences_expert(destination_file)  # id_pubmed : [sentence1, sentence2, etc.]

    for key, items in dict_id_sentence.items():
        root = ET.Element('document', id=key)
        sentence_number = 0

        for item in items:
            entity_number = 0
            #entities_sentence = get_entities_sentence(destination_file, item)
            entities_sentence = get_entities_sentence_expert(destination_file, item)  # for expert

            doc = ET.SubElement(root, 'sentence', id=key + '.s' + str(sentence_number), text=item)
            save_entities = []
            save_pairs = []

            for pair in entities_sentence:  # sentence : [[entity1, id, char1, char2], [entity2, id, char3, char4],
                # relation], [[entity1, id, char1, char2], [entity2, id, char3, char4],
                # relation], etc.]

                pair_list = []

                entity_1 = pair[0]
                entity_2 = pair[1]

                if entity_1 not in save_entities:
                    save_entities.append(entity_1)

                if entity_2 not in save_entities:
                    save_entities.append(entity_2)

                pair_list.append(entity_1)
                pair_list.append(entity_2)

                pair_list_sorted = sorted(pair_list, key=lambda x: int(x[2]))
                pair_list_sorted.append(pair[2])

                save_pairs.append(pair_list_sorted)

            save_entities_sorted = sorted(save_entities, key=lambda x: int(x[2]))
            associated_entity_number = {}

            for entity in save_entities_sorted:

                if entity[1].startswith('HP'):

                    ET.SubElement(doc, 'entity', id=key + '.s' + str(sentence_number) + '.e' + str(entity_number),
                                  charOffset=entity[2] + '-' + entity[3], type='HP', text=entity[0],
                                  ontology_id=entity[1])

                else:

                    ET.SubElement(doc, 'entity', id=key + '.s' + str(sentence_number) + '.e' + str(entity_number),
                                  charOffset=entity[2] + '-' + entity[3], type='GENE', text=entity[0],
                                  ontology_id=entity[1])

                associated_entity_number[repr(entity)] = key + '.s' + str(sentence_number) + '.e' + str(entity_number)
                entity_number += 1

            # save_pairs = [[['PDYN', '5173', '64', '68'], ['epilepsy', 'HP_0001250', '111', '119'], 'False'], [['PDYN', '5173', '32', '68'], ['epilepsy', 'HP_0001250', '111', '119'], 'False']]

            pairs_sorted = sorted(save_pairs, key=lambda x: int(x[0][2]))
            pair_number = 0

            for pair in pairs_sorted:

                ET.SubElement(doc, 'pair', id=key + '.s' + str(sentence_number) + '.p' + str(pair_number),
                              e1=associated_entity_number[repr(pair[0])],
                              e2=associated_entity_number[repr(pair[1])],
                              relation=pair[2].lower())

                pair_number += 1
            sentence_number += 1

        output_file = open(destination_path + key + '.xml', 'w', encoding='utf-8')
        # output_file.write('<?xml version="1.0" encoding="UTF-8"?>')

        #root = xml_file_go(dict_gene_go_name, dict_gene_go_id, root)
        output_file.write(prettify(root))

        output_file.close()

    return


#get_original_test_set_xml('data/original_test_set.tsv', 'data/original_test_set_consensus.tsv', 'corpora/original_test/', 'data/gene2go')
#get_original_test_set_xml('data/original_test_set.tsv', 'data/original_dataset_70.tsv', 'corpora/original_train/', 'data/gene2go')
#get_original_test_set_xml('data/original_test_set.tsv', 'data/expert.tsv', 'corpora/expert_test/', 'data/gene2go')

def all_in_one(base_dir, destination_file, destination_path, nomenclature='relation'):
    """

    :param nomenclature:
    :param base_dir:
    :param destination_file:
    :param destination_path:
    :return:
    """

    base_root = ET.Element('document', id=destination_file)
    sentence_number = 0

    for f in os.listdir(base_dir):

        reader = open(base_dir + '/' + f, 'r', encoding='utf-8')
        content = reader.read().replace('</sup>', 'AAAA').replace('<sup>', 'AAA').replace('</b>', 'AAAA')\
            .replace('<b>', 'AAA').replace('</i>', 'AAAA').replace('<i>', 'AAA').replace('</sub>', 'AAAA')\
            .replace('<sub>', 'AAA')

        root = ET.fromstring(content)

        #tree = ET.parse(base_dir + '/' + f)
        #root = tree.getroot()

        for sentence in root:
            sentence_text = sentence.get('text')
            all_pairs = sentence.findall('pair')

            for pair in all_pairs:
                doc = ET.SubElement(base_root, 'sentence', id='s' + str(sentence_number), text=str(sentence_text.encode('ascii', 'ignore'))[2:-1])

                id_1 = pair.get('e1')
                id_2 = pair.get('e2')

                for e in sentence.findall('entity'):

                    if e.get('id') == id_1:
                        if e.get('type') == 'GO':
                            e_type = 'GO'
                        else:
                            e_type = 'HP'

                        ET.SubElement(doc, 'entity', id='s' + str(sentence_number) + '.e0',
                                      charOffset=e.get('charOffset'), type=e_type, text=e.get('text'), ontology_id=e.get('ontology_id'))

                    elif e.get('id') == id_2:
                        if e.get('type') == 'GO':
                            e_type = 'GO'
                        else:
                            e_type = 'HP'

                        ET.SubElement(doc, 'entity', id='s' + str(sentence_number) + '.e1',
                                      charOffset=e.get('charOffset'), type=e_type, text=e.get('text'), ontology_id=e.get('ontology_id'))

                ET.SubElement(doc, 'pair', id='s' + str(sentence_number) + '.p0', e1='s' + str(sentence_number) + '.e0',
                              e2='s' + str(sentence_number) + '.e1', relation=pair.get(nomenclature))

                sentence_number += 1

    all_file = open(destination_path + '/' + destination_file + '.xml', 'w', encoding='utf-8')
    all_file.write(prettify(base_root))
    all_file.close()

    return


#all_in_one('corpora/consensus_test/', 'consensus_test', 'corpora/')
#all_in_one('corpora/original_train/', 'original_train', 'corpora/')
#all_in_one('corpora/original_test/', 'original_test', 'corpora/')
#all_in_one('corpora/expert_test/', 'expert_test', 'corpora/')

# RUN ####

# def main():
#     """Creates an xml file for each abstract
#     """
#
#     #xml_file('data/original_dataset_70.tsv', 'data/batch_results_70.csv', 'data/validation_set.csv', 'data/gene2go', 'corpora/amazon_train/')
#     #xml_file('data/original_dataset_30.tsv', 'data/batch_results_30_consensus.csv', 'data/validation_set.csv', 'data/gene2go', 'corpora/consensus_test/', test=True)
#
#     return
#
#
# # python3 src/parser_csv.py
# if __name__ == "__main__":
#     main()
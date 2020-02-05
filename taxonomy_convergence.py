import pandas as pd



def get_types(taxon_ranks,input_data,field_name):
    replicons = {}

    num_ranks = len(taxon_ranks)

    for i,row in input_data.T.iteritems():
        if row[field_name] == '-':
            continue
        replicons_types = str(row[field_name]).split(',')

        for k in range(0,num_ranks):
            taxon = taxon_ranks[k]
            taxon_name = str(row[taxon])
            if taxon_name != 'nan':
                names = [taxon_name]
            j = k

            while taxon_name == 'nan':
                names.append("{}".format(taxon_ranks[j]))
                j-=1
                taxon_name = row[taxon_ranks[j]]


            if len(names) > 1:
                taxon_name = "__incerta_sedis__".join(names)


            for rep in replicons_types:

                if rep not in replicons:
                    replicons[rep] = {}

                if not taxon  in replicons[rep]:
                    replicons[rep][taxon] = {}

                if not taxon_name in replicons[rep][taxon]:
                    replicons[rep][taxon][taxon_name] = 0


                replicons[rep][taxon][taxon_name]+=1
    return replicons



def process_results(taxon_ranks,replicons):
    num_ranks = len(taxon_ranks)
    ranges = {}
    for i in range(0,num_ranks):
        taxon = taxon_ranks[i]

    for rep in replicons:
        convergence_found = False
        for i in range(0,num_ranks):
            taxon = taxon_ranks[i]
            if len(replicons[rep][taxon]) > 1:

                convergence_found = True

                j = i - 1
                convergence_rank = taxon_ranks[j]
                convergence_name = list(replicons[rep][taxon_ranks[j]].keys())[0]
                convergence_count = replicons[rep][taxon_ranks[j]][convergence_name]

                while 'incerta_sedis' in convergence_name:
                    j-=1
                    convergence_rank = taxon_ranks[j]
                    convergence_name = list(replicons[rep][taxon_ranks[j]].keys())[0]
                    convergence_count = replicons[rep][taxon_ranks[j]][convergence_name]

                break
        if not convergence_found:
            convergence_rank = 'species'
            convergence_name = list(replicons[rep]['species'].keys())[0]
            convergence_count = replicons[rep]['species'][convergence_name]

        ranges[rep] = {
            'convergence_rank': convergence_rank,
            'convergence_name': convergence_name,
            'convergence_count': convergence_count
        }


    return ranges

def determine_overall_convergence(taxon_ranks,input_data,replicon_types,relaxase_acs,cluster_codes):

    num_ranks = len(taxon_ranks)



    for i,row in input_data.T.iteritems():
        overall_ranks = {}
        accession = row['accession']
        record_replicons_types = str(row['rep_type(s)']).split(',')
        record_relaxase_acs = str(row['relaxase_type_accession(s)']).split(',')

        record_cluster_code = str(row['mash_neighbor_cluster'])
        convergence_rank_rep = '-'
        convergence_name_rep = '-'
        top_rank_index = 99999
        if (row['rep_type(s)'] != '-'):
            for rep in record_replicons_types:
                convergence_rank_rep = replicon_types[rep]['convergence_rank']
                convergence_name_rep = replicon_types[rep]['convergence_name']
                convergence_rank_index = taxon_ranks.index(convergence_rank_rep)
                if not convergence_rank_rep in overall_ranks:
                    overall_ranks[convergence_rank_rep] = {}
                overall_ranks[convergence_rank_rep ][convergence_name_rep] = ''

                if convergence_rank_index < top_rank_index:
                    top_rank = convergence_rank_rep
                    top_taxa = convergence_name_rep
                    top_rank_index = convergence_rank_index

            convergence_rank_rep = top_rank
            convergence_name_rep = top_taxa

        convergence_rank_rel = '-'
        convergence_name_rel = '-'
        top_rank_index = 99999
        if (row['relaxase_type_accession(s)'] != '-'):
            for rel in record_relaxase_acs:
                convergence_rank_rel = relaxase_acs[rel]['convergence_rank']
                convergence_name_rel = relaxase_acs[rel]['convergence_name']
                convergence_rank_index = taxon_ranks.index(convergence_rank_rel)
                if not convergence_rank_rel in overall_ranks:
                    overall_ranks[convergence_rank_rel] = {}
                overall_ranks[convergence_rank_rel][convergence_name_rel] = ''

                if convergence_rank_index < top_rank_index:
                    top_rank = convergence_rank_rel
                    top_taxa = convergence_name_rel
                    top_rank_index = convergence_rank_index
            convergence_rank_rel = top_rank
            convergence_name_rel = top_taxa

        convergence_rank_code = cluster_codes[record_cluster_code]['convergence_rank']

        if not convergence_rank_code in overall_ranks:
            overall_ranks[convergence_rank_code] = {}
        convergence_name_code = cluster_codes[record_cluster_code]['convergence_name']


        overall_ranks[convergence_rank_code][convergence_name_code ] = ''
        #print(overall_ranks)

        overall_convergence_rank = 'species'
        overall_convergence_name = row['species']

        for i in range(0,num_ranks):
            rank = taxon_ranks[i]
            if rank not in overall_ranks:
                continue

            overall_convergence_name = row[rank]
            overall_convergence_rank = rank
            break


        print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(accession,overall_convergence_rank,overall_convergence_name,convergence_rank_rep,convergence_name_rep,convergence_rank_rel,convergence_name_rel,convergence_rank_code,convergence_name_code))




def main():
    taxon_ranks = ['superkingdom',
    'phylum',
    'subphylum',
    'superclass',
    'class',
    'subclass',
    'infraclass',
    'superorder',
    'order',
    'suborder',
    'infraorder',
    'superfamily',
    'family',
    'subfamily',
    'genus',
    'species']


    input_data = pd.read_csv('/Users/jrobertson/PycharmProjects/__mobcluster/__2019-12-results/Taxonomy_NCBI_Plasmids.txt', header=0,sep="\t")

    replicon_types = process_results(taxon_ranks,get_types(taxon_ranks, input_data, 'rep_type(s)'))
    replicon_acs = process_results(taxon_ranks,get_types(taxon_ranks, input_data, 'rep_type_accession(s)'))
    relaxase_types = process_results(taxon_ranks, get_types(taxon_ranks, input_data, 'relaxase_type(s)'))
    relaxase_acs = process_results(taxon_ranks, get_types(taxon_ranks, input_data, 'relaxase_type_accession(s)'))
    cluster_codes = process_results(taxon_ranks, get_types(taxon_ranks, input_data, 'mash_neighbor_cluster'))
    determine_overall_convergence(taxon_ranks, input_data, replicon_types, relaxase_acs, cluster_codes)




main()
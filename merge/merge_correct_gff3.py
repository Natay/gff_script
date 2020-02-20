# This script merges stringtie annotated with augustus annoatated gff3.
# A transcript that overlaps with more than one gene is removed.
# usage: python merge_correct_gff3.py <aug_ann.gff3> <stringtie_uniq_ann.gff3> <species>

import sys, csv, copy
from collections import defaultdict


def extract_attr(attr_str, attr_name):
    attrs = attr_str.split(";")
    feat_id = [val for val in attrs if val.startswith(attr_name)][0]
    feat_id = feat_id.replace(attr_name, "")
    return feat_id


def extract_tid(attr):
    tid = attr.split(";")[1].split(" ")
    tid = list(filter(lambda x: x != '', tid))[1]
    tid = tid.replace("\"", "")
    return tid


def read_gff(gff, feat="gene"):
    """
    Reads augustus annotated gff3 file.
    """
    ann_store = defaultdict(list)
    stream = csv.reader(open(gff), delimiter="\t")

    for row in stream:

        if row[0].startswith("#"):
            continue
        if row[2] == feat:
            chrom, start, end, strand, attr = row[0], row[3], row[4], row[6], row[8]
            feat_id = extract_attr(attr, "ID=")
            key = ":".join([chrom, start, end, strand, feat_id])
            ann_store[key].append(row)
            continue

        ann_store[key].append(row)
        print(key)
        1/0
    return ann_store


def read_ann(fname):
    """
    Stores augustus annotated genes
    """
    ann_store = defaultdict(list)
    stream = csv.DictReader(open(fname), delimiter="\t")
    for row in stream:
        key = ":".join([row['chrom'],row['start'],row['end'], row['strand']])
        ann_store[key] = row
    return ann_store


def store_recp_map(fname, species):
    """
    Reads a tab-delimited idmap file and returns a dictionary of list.
    Eg: giraffe_recp_protein_map_stringtie.txt
    """
    trans_store = defaultdict(dict)
    stream = csv.DictReader(open(fname), delimiter="\t")
    for row in stream:
        protein_id = species + "_pid"
        pid = row.get(protein_id)
        trans_store[pid] = row
    return trans_store


def is_overlap(chrom1, start1, end1, strand1, chrom2, start2, end2, strand2):
    """
    Checks coordinate overlap.
    """
    if chrom1 != chrom2:
        return False
    #if strand1 != strand2:
    #    return False

    cov = 0
    len1 = abs(end1 - start1) + 1
    len2 = abs(end2 - start2) + 1

    # check overlap - within a transcript
    if start1 >= start2 and end1 <= end2:
        cov = ((end1 - start1) / len1) * 100

    if start1 <= start2 and end1 >= end2:

        cov = ((end2 - start2) / len2) * 100

    if start1 < start2 and (end1 >= start2 and end1 < end2):

        # overlap left
        cov = ((end1 - start2) / len1) * 100
    if start1 > start2 and (start1 < end2 and end1 > end2):
        # overlap right
        cov = ((end2 - start1) / len1) * 100
    if chrom1 == chrom2 and cov >= 60:
        return True
    return False


def check_gene_overlap(store, transcript):
    """
    Check if a stringtie transcript overlap with an augustus predicted gene.
    """
    gene_overlap = list()

    tr_chrom, tr_start, tr_end, tr_strand, tr_id = transcript.split(":")
    tr_start, tr_end = int(tr_start), int(tr_end)

    for gene in store.keys():
        gene_chrom, gene_start, gene_end, gene_strand, gene_id = gene.split(":")
        gene_start, gene_end = int(gene_start), int(gene_end)

        if is_overlap(tr_chrom, tr_start, tr_end, tr_strand, gene_chrom, gene_start, gene_end, gene_strand):
            gene_overlap.append(gene)
    return gene_overlap


def get_attributes(curr_attr):
    uid = extract_attr(curr_attr, "ID=")
    gid = ".".join(uid.split(".")[:-1])
    gid = gid if gid else extract_attr(curr_attr, "gene_id=")
    gname = extract_attr(curr_attr, "gene_name=")
    best_match = extract_attr(curr_attr, "best_match=")
    source = best_match
    evidence = get_evidence(source)
    pident, qcovs = get_evidence_vals(evidence, curr_attr)

    return uid, gid, gname, best_match, source, evidence, pident, qcovs


def modify_stringtie_gene_attr(gid, gene_id, gname, source, qcovs, pident,evidence):

    if evidence == "stringtie" or evidence == "augustus":
        new_attr = f'ID={gid};gene_id={gid};gene_name={gene_id};Name={gene_id};best_match=stringtie'
    else:
        new_attr = f'ID={gid};gene_id={gid};gene_name={gname};Name={gname};best_match={source};query_coverage={qcovs};percentage_identity={pident}'

    return new_attr


def modify_gene_attr(gid, gname, source, qcovs, pident,evidence):
    if evidence == "stringtie" or evidence == "augustus":
        new_attr = f'ID={gid};gene_id={gid};gene_name={gid};Name={gid};best_match=stringtie'
    else:
        new_attr = f'ID={gid};gene_id={gid};gene_name={gname};Name={gname};best_match={source};query_coverage={qcovs};percentage_identity={pident}'

    return new_attr


def modify_transcript_attr(uid, tid, tidx, gid, gname, source, qcovs, pident, evidence, term, seen_tuids):
    parent = gid
    trans_no = tid.split(".")[1]
    uid_mod = gid + ".t" + str(tidx)

    while uid_mod in seen_tuids:
        tidx+=1
        uid_mod = gid + ".t" + str(tidx)

    if evidence == "stringtie":
        trans_id = gid + ".t" + str(tidx) if term != "gene"  else uid_mod
        new_attr = f'ID={uid_mod};Parent={parent};original_id={uid};transcript_id={trans_id};gene_name={parent};Name={trans_id};best_match={evidence}'
    elif evidence == "locus":
        trans_id = gname + ".t" + str(tidx) if term != "gene"  else ".".join([gname, trans_no])
        new_attr = f'ID={uid_mod};Parent={parent};original_id={uid};transcript_id={trans_id};gene_name={gname};Name={trans_id};best_match={evidence}'

    else:
        trans_id = gname + ".t" + str(tidx) if term != "gene"  else ".".join([gname, trans_no])
        if term == "stringtie_transcript"  or term == "parent_change":
            new_attr = f'ID={uid_mod};Parent={parent};original_id={uid};transcript_id={trans_id};gene_name={gname};Name={trans_id};best_match={source};query_coverage={qcovs};percentage_identity={pident}'
        else:
            new_attr = f'ID={uid};Parent={parent};original_id={uid};transcript_id={trans_id};gene_name={gname};Name={trans_id};best_match={source};query_coverage={qcovs};percentage_identity={pident}'
    return new_attr


def modify_attr(attr, tid, tidx, gname, term):
    parent = tid
    trans_no = tid.split(".")[1]

    if gname is None or gname == "":
        new_attr = f'Parent={parent};transcript_id={parent}'
    else:
        trans_id = gname + ".t" + str(tidx) if term != "gene"  else ".".join([gname, trans_no])
        new_attr = f'Parent={parent};transcript_id={trans_id};gene_name={gname}'

    if attr.startswith("ID"):
        feat = attr.split(";")[0].split(".")[-1]
        id_ = "ID=" + ".".join([parent, feat])
        new_attr = ";".join([id_, new_attr])

    return new_attr


def annotate_existing(attr, count, existing_vals, ann_term, parent_gene, seen_ids=set()):
    vals1 = copy.deepcopy(existing_vals)

    #uid, gid, gname, best_match, source, evidence, pident, qcovs = get_attributes(attr)
    d = dict()
    #seen_ids = set()

    for item in vals1:

        if item[2] == "gene":
            gid = extract_attr(item[8], "ID=")  # keeping the original ID
            gene_attr = modify_gene_attr(gid, gname, source, qcovs, pident, evidence)
            item[8] = gene_attr

        elif item[2] == "transcript":

            tuid, tgid, tgname, tbest_match, tsource, tevidence, tpident, tqcovs = get_attributes(item[8])
            d[tuid] = parent_gene

            if ann_term == "transcript":
                uid, gid, gname, best_match, source, evidence, pident, qcovs = tuid, gid, gname, tbest_match, tsource, tevidence, tpident, tqcovs
                tid = uid

            if ann_term == "gene":
                uid, gid, gname, best_match, source, evidence, pident, qcovs = tuid, tgid, gname, tbest_match, tsource, tevidence, tpident, tqcovs
                tid = uid

            if ann_term == "parent_change":
                gid = parent_gene
                tid = uid

            if ann_term == "transcript_locus":
                source, evidence = "locus", "locus"
                gid = parent_gene
                tid = tuid

            tran_attr = modify_transcript_attr(uid, tid, count, gid, gname, source, qcovs, pident, evidence,
                                               ann_term, seen_ids)
            item[8] = tran_attr
            d[tuid] = extract_attr(tran_attr, "ID=")
            seen_ids.add(d[tuid])

        else:
            #TODO: what is the point of d
            pid = d[tid]
            new_attr = modify_attr(item[8], pid, count, gname, ann_term)
            item[8] = new_attr

    return vals1


def get_evidence(source):
    if source not in ["stringtie", "locus", "augustus"]:
        evidence = "species"
    else:
        evidence = source
    return evidence


def get_evidence_vals(evidence, attrs):
    """
    returns pident, qcovs
    """
    pident, qcovs = 0, 0
    if evidence == "species":
        pident = extract_attr(attrs, "percentage_identity=")
        qcovs = extract_attr(attrs, "query_coverage=")
    return pident, qcovs


def get_merged_keys(merged):
    store = list()
    for item in merged:
        if item[2] != "gene":
            continue
        chrom, start, end, strand = item[0], item[3], item[4], item[6]
        key = ":".join([chrom, start,end,strand])
        store.append(key)
    return store


def get_keys(vals):
    key_vals = [":".join(v.split(":")[:-1]) for v in vals]
    return key_vals


def get_unchaged_keys(aug_store, current):
    """
    Extract augustus_store keys that are not in current set of keys.
    """

    prefixes = [":".join(key.split(":")[:-1]) for key in aug_store.keys()]
    stable =[k for k in prefixes if k not in current]
    unchanged = [k for k in aug_store.keys() if ":".join(k.split(":")[:-1]) in stable ]

    return unchanged


def modify_transcript_parent(genes):
    """
    Make transcript parent = gene uid
    """
    vals  = copy.deepcopy(genes)
    for gene, transcripts in vals.items():
        for t in transcripts:
            if t[2] not in ["gene", "transcript"]:
                continue
            if t[2] == "gene":
                guid = extract_attr(t[8], "ID=")
                continue
            if t[2] == "transcript":
                parent = "Parent=" + extract_attr(t[8], "Parent=")
                mod = t[8].replace(parent, f"Parent={guid}")
                t[8] = mod
                continue
    return vals


def modify_stringtie_transcript_count(transcripts, seen_tids):
    trans = copy.deepcopy(transcripts)
    tgname = ""
    #seen_trans = set()
    parent_dict = dict()
    for gene, values in trans.items():
        count = 1
        for item in values:
            if item[2] == "gene":
                continue

            elif item[2] == "transcript":

                tuid, tgid, tgname, tbest_match, tsource, tevidence, tpident, tqcovs = get_attributes(item[8])
                parent_gene = extract_attr(item[8], "Parent=") # modify transcript parent
                tattrs= modify_transcript_attr(tuid, tuid, count, parent_gene, tgname, tsource, tqcovs, tpident, tevidence,
                                               "stringtie_transcript", seen_tids )
                item[8] = tattrs
                count+=1
                parent_dict[tuid] = extract_attr(tattrs, "ID=")
            else:

                tuid = extract_attr(item[8], "Parent=")
                feat_parent = parent_dict[tuid]
                new_attr = modify_attr(item[8], feat_parent, count, tgname, "stringtie_transcript")
                item[8] = new_attr

    return trans


def get_transcript_count(elms):
    """
    Extract the no. of transcripts.
    """
    tcount= 0
    for elm  in elms:
        if elm[2] == "transcript":
            tcount+=1
    return tcount


def create_merged_gff(aug_store, transcript_store):
    new_trans_store = list()
    new_gene_store = defaultdict(list)
    merged_store = list()
    seen, seen_tids = set(), set()

    i, tx = 0, 0

    prev_gene = ""

    for transcript_key, transcript_vals in transcript_store.items():

        trans_attrs = transcript_vals[0][8]
        tr_best_match = extract_attr(trans_attrs, "best_match=")

        #  check transcript overlap with augustus annotated gene
        overlap_store = check_gene_overlap(aug_store, transcript_key)

        # CASE 1 : # transcript overlaps with more than one gene - remove transcript.
        if len(overlap_store) > 1:
            # del transcript_store_corrected[transcript_key]
            continue

        # Overlap with an augustus gene
        if len(overlap_store) == 1:
            gene_key = overlap_store[0]

            # Discard if the overlapped gene and transcript do not have the same orientation
            gene_strand = gene_key.split(":")[3]
            transcript_strand = transcript_key.split(":")[3]

            if gene_strand != transcript_strand:
                continue

            gene_attr = aug_store[gene_key][0][8]
            gene_best_match = extract_attr(gene_attr, "best_match=")
            tcount = get_transcript_count(aug_store[gene_key])
            tx = tx+1 if gene_key == prev_gene else 1
            tcount +=tx
            prev_gene = gene_key

            # CASE2 : transcript is annotated and gene is annotated -
            # Nothing to do except correct transcript suffix numbers and assigning correct parent id
            if tr_best_match != "stringtie" and gene_best_match != "augustus":
                #print("CASE2")
                annotate = "transcript"  # keeping the annotations from transcript itself, change the  parent to augustus gene
                parent_gene = extract_attr(gene_attr,"ID=")

                modified_transcript = annotate_existing(gene_attr, tcount, transcript_vals, annotate, parent_gene,
                                                        seen_tids)
                if gene_key not in seen:
                    merged_store.extend(aug_store[gene_key])
                    seen.add(gene_key)
                merged_store.extend(modified_transcript)
                continue

            # CASE3 : transcript is unannotated and gene is unannotated
            # Nothing to do except correct transcript suffix numbers and assign correct parent id
            if tr_best_match == "stringtie" and gene_best_match == "augustus":
                #print("CASE3")
                annotate = "transcript_locus" #
                parent_gene = extract_attr(gene_attr,"ID=")
                modified_transcript = annotate_existing(gene_attr, tcount, transcript_vals, annotate, parent_gene,
                                                        seen_tids)
                if gene_key not in seen:
                    merged_store.extend(aug_store[gene_key])
                    seen.add(gene_key)
                merged_store.extend(modified_transcript)
                continue

            # CASE 4 : transcript is annotated gene is unannotated - annotate gene with transcript ann
            # chaange transcript uid and parent
            if tr_best_match != "stringtie"  and gene_best_match == "augustus":
                #print("CASE4")
                annotate = "gene"
                gene_vals = aug_store[gene_key]
                parent_gene = extract_attr(gene_vals[0][8], "ID=")
                modified_gene  =  annotate_existing(trans_attrs, tcount, gene_vals, annotate, parent_gene,seen)
                modified_transcript = annotate_existing(trans_attrs, tcount, transcript_vals, "parent_change",
                                                        parent_gene,seen_tids)

                if gene_key not in seen:
                    merged_store.extend(modified_gene)
                    seen.add(gene_key)
                merged_store.extend(modified_transcript)
                continue

            # CASE 5 : transcript is unannotated and gene is annotated - annotate transcript with gene ann
            if tr_best_match == "stringtie" and gene_best_match != "augustus":
                #print("CASE5")
                annotate = "transcript_locus"
                parent_gene = extract_attr(gene_attr,"ID=")
                modified_transcript = annotate_existing(gene_attr, tcount, transcript_vals, annotate,parent_gene,
                                                         seen_tids)
                if gene_key not in seen:
                    merged_store.extend(aug_store[gene_key])
                    seen.add(gene_key)
                merged_store.extend(modified_transcript)
                continue

        # No overlap - add gene row to the largest transcript
        if len(overlap_store) == 0:
            #print("NO overlap")
            # # check if it is overlaps with any other transcript seen so far.
            # # specify gene length as the length of the largest transcript.
            tr_chrom, tr_start, tr_end, tr_strand, tr_id = transcript_key.split(":")
            tr_start, tr_end = int(tr_start), int(tr_end)
            gchrom, gstart, gend, gstrand = tr_chrom, tr_start, tr_end, tr_strand

            if new_trans_store:
                for item in new_trans_store:

                    chrom, start, end, strand, tid = item
                    start, end = int(start), int(end)

                    if is_overlap(chrom, start, end, strand, tr_chrom, tr_start, tr_end, tr_strand):
                        gstart = tr_start if tr_start < start else start
                        gend = tr_end if tr_end > end else end

            best_match =extract_attr(trans_attrs, "best_match=")
            source = best_match
            tid = extract_attr(trans_attrs, "ID=")
            gid = ".".join(tid.split(".")[:-1])
            gname = extract_attr(trans_attrs, "gene_name=")
            evidence = get_evidence(best_match)
            if evidence == "augustus" or evidence == "stringtie":
                gname = gid
            pident, qcovs = get_evidence_vals(evidence, trans_attrs)

            #  There are transcripts which stringtie thinks should come from the same gene.
            # But augustus predicted two genes in that location instead of one which also got annotated.
            # So, assign uniq ids to these genes.
            stringtie_gene_id =gid
            i+=1
            gid = gid +".gene" + str(i)
            gene_attr = modify_stringtie_gene_attr(gid, stringtie_gene_id, gname, source, qcovs, pident, evidence)
            #gene_attr = modify_gene_attr(gid, gname, source, qcovs, pident, evidence)

            gene_row = [gchrom, "StringTie", "gene" , gstart, gend, ".", gstrand, ".", gene_attr]

            if gname not in new_gene_store:
                new_gene_store[gname].append(gene_row)
            else:
                new_gene_store[gname][0] = gene_row

            new_gene_store[gname].extend(transcript_vals)
            new_trans_store.append((gchrom, gstart, gend, gstrand, transcript_key))

            continue

    # change the parent attribute of the strigtie-specific transcripts
    # to the newly assigned gene-id.
    new_gene_store = modify_transcript_parent(new_gene_store)
    new_gene_store = modify_stringtie_transcript_count(new_gene_store, seen_tids)
    return merged_store, new_gene_store


if __name__ == "__main__":
    aug_gff = "test_aug.gff3"
    stringtie_gff = "stringtie_uniq_ann.gff3"
    species = "giraffe"

    aug_gff_store = read_gff(aug_gff, "gene")
    transcript_gff_store = read_gff(stringtie_gff, "transcript")

    merged, added = create_merged_gff(aug_gff_store, transcript_gff_store)

    current_keys = list()
    if merged:
        merged_keys = get_merged_keys(merged)
        current_keys.extend(merged_keys)
    if added:
        added_keys =get_keys(added.keys())
        current_keys.extend(added_keys)

    if not current_keys:
        print ("Nothing added from stringtie file")
        current_keys = aug_gff_store.keys()
        # print augustus and exit
        for k , v in aug_gff_store.items():
            for item in v:
                item = list(map(str, item))
                print("\t".join(item))
        sys.exit()

    # Get augustus genes that did not have any intersection with stringtie, ie, augustus uniq
    unchanged = get_unchaged_keys(aug_gff_store, current_keys)

    #print header
    print("##gff-version 3")

    if unchanged:
        # print unchanged
        for key in unchanged:
            for item in aug_gff_store[key]:
                item = list(map(str, item))
                print("\t".join(item))

    if merged:
        # print merged:
        for item in merged:
            item = list(map(str, item))
            print("\t".join(item))

    if added:
        # print additional from stringtie
        for gene, childeren in added.items():
            for child in childeren:
                child= list(map(str, child))
                print("\t".join(child))

import csv
import copy

TEST_AUGUSTUS = "aug.gff3"
TEST_STRINGTIE = "stringtie.gff3"


def parse_row(row):
    # chrom, start, end, strand, attr, feat
    return row[0], row[3], row[4], row[6], row[8], row[2]


def is_augustus(data):
    return data['t_best'] != "stringtie" and data['g_best'] != "augustus"


def dict_to_str(data, delim=";"):

    lst = [f"{k}={v}" for k, v in data.items()]
    stri = f"{delim}".join(lst)

    return stri


def parse_attrs(attr_str, attr_name):
    # Find start and ending index of the attribute we want
    if not attr_str:
        return ""
    start_idx = attr_str.find(attr_name)
    end_idx = attr_str[start_idx:].find(";")

    # Slice string to the attribute we want
    value = attr_str[start_idx:end_idx]
    value = value.replace(f"{attr_name}=", "")

    return value


def parse_file(gff, target="gene"):
    ann_store = dict()
    meta_data = dict()
    stream = csv.reader(open(gff), delimiter="\t")
    key = None
    for row in stream:
        if row.startswith("#"):
            continue

        chrom, start, end, strand, attr, feat = parse_row(row)
        if feat == target:
            key = parse_attrs(attr, "ID")
            # Fill meta data one for each key
            meta_data[key] = dict(chrom=chrom, start=start, end=end, attr=attr, strand=strand)
        if key:
            ann_store.setdefault(key, []).append(row)

    return ann_store, meta_data


def calculate_coverage(tmeta, gmeta):

    ch1, ch2 = tmeta['chrom'], gmeta['chrom']
    en1, en2 = tmeta['end'], gmeta['end']
    st1, st2 = tmeta['start'], gmeta['start']

    cov = 0
    len1 = abs(en1 - st1) + 1
    len2 = abs(en2 - st2) + 1

    within = st1 >= st2 and en1 <= en2
    within2 = st1 <= st2 and en1 >= en2
    # Parameters to evaluate left and right transcript coverage.
    left = st1 < st2 and en1 >= st2 and en1 < en2
    right = st1 > st2 and st2 < en2 and en1 > en2

    if within:
        cov = ((en1 - st1) / len1) * 100
    if within2:
        cov = ((en2 - st2) / len2) * 100
    if left:
        cov = ((en1 - st2) / len1) * 100
    if right:
        cov = ((en2 - st1) / len1) * 100

    return cov, ch1, ch2


def is_overlap(tmeta, gmeta):
    # Parse relevant coordinates from meta data
    coverage, ch1, ch2 = calculate_coverage(tmeta=tmeta, gmeta=gmeta)

    if ch1 != ch2:
        return False
        # if strand1 != strand2:
        #    return False

    if ch1 == ch2 and coverage >= 60:
        return True

    return False


def check_overlap(gmeta, tmeta):
    overlap_store = []
    for gene, gene_meta in gmeta.items():
        # Append gene to list of overlapping genes
        if is_overlap(tmeta=tmeta, gmeta=gmeta):
            overlap_store.append(gene)
    return overlap_store


def get_count(rows, target_str, target_idx):

    # Filter for rows that have the target_str in the index ( target_idx )
    target_lst = list(filter(lambda l: l[target_idx] == target_str, rows))
    count = len(target_lst)
    return count


def parse_data(gmeta_data, trans_meta, gstore):
    # Keep track of genes overlapping with this transcript.

    trans_strand = trans_meta['strand']

    overlap_store = check_overlap(gmeta_data, trans_meta)
    is_gene_overlap = len(overlap_store) == 1
    is_overlapped = len(overlap_store) > 1
    no_overlap = len(overlap_store) == 0
    g_key = overlap_store[0] if is_gene_overlap else ""
    gene_meta = gmeta_data.get(g_key, {})
    gene_strand = gene_meta.get('strand', '')

    t_best = parse_attrs(trans_meta['attr'], "best_match")
    g_best = parse_attrs(gene_meta.get('attr', ""), "best_match")

    data = dict(trans_strand=trans_strand, overlap_store=overlap_store,
                is_overlapped=is_overlapped,
                no_overlap=no_overlap,
                is_gene_overlap=is_gene_overlap,
                g_key=g_key, gene_meta=gene_meta,
                gene_strand=gene_strand, gene=gstore[g_key], t_best=t_best, g_best=g_best)
    return data


def get_evidence(source):
    if source in ["stringtie", "locus", "augustus"]:
        return source
    evidence = "species"
    return evidence


def get_evidence_vals(evidence, attrs):
    """
    returns pident, qcovs
    """
    pident, qcovs = 0, 0
    if evidence == "species":
        pident = parse_attrs(attrs, "percentage_identity")
        qcovs = parse_attrs(attrs, "query_coverage")
    return pident, qcovs


def modify_gene_attrs(gid, attrs_dict):
    is_augus = attrs_dict['is_augus']
    name = attrs_dict['gname']
    source = attrs_dict['source']
    qcovs = attrs_dict['qcovs']
    pident = attrs_dict['pident']
    modified = dict(ID=gid, gene_id=gid, gene_name=gid, Name=gid,
                    best_match=source)
    if not is_augus:
        update_dict = dict(gene_name=name, Name=name, best_match=source,
                           query_coverage=qcovs, percentage_identity=pident)
        modified.update(update_dict)

    modified_str = dict_to_str(data=modified)

    return modified_str


def modify_transcript_attr(uid, tid, tidx, gid, gname, source, qcovs, pident, evidence, term, seen_tuids):
    parent = gid
    trans_no = tid.split(".")[1]
    uid_mod = gid + ".t" + str(tidx)

    while uid_mod in seen_tuids:
        tidx += 1
        uid_mod = gid + ".t" + str(tidx)

    trans_id = gid + ".t" + str(tidx)
    attrs_dict = dict()
    if evidence == "stringtie":
        trans_id = gid + ".t" + str(tidx) if term != "gene" else uid_mod
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


def get_attrs_dict(full_attr):
    uid = parse_attrs(full_attr, "ID")
    gid = ".".join(uid.split(".")[:-1])
    gid = gid if gid else parse_attrs(full_attr, "gene_id")
    gname = parse_attrs(full_attr, "gene_name")
    source = parse_attrs(full_attr, "best_match")
    evidence = get_evidence(source)
    pident, qcovs = get_evidence_vals(evidence, full_attr)
    # Is augustus
    is_augus = evidence == "stringtie" or 'evidence' == "augustus"

    data = dict(uid=uid, gid=gid, gname=gname, best_match=source, source=source,is_augus=is_augus,
                evidence=evidence, pident=pident, qcovs=qcovs)
    return data


def annotate_transcript(parent_gene, term, gene_attrs_dict, trans_attr, seen_ids=set(), collect=dict()):
    trans_attrs = get_attrs_dict(trans_attr)
    tuid = trans_attrs["uid"]
    # Add transcript to collect dict.
    collect[tuid] = parent_gene

    uid, gid = gene_attrs_dict["uid"], gene_attrs_dict["gid"]
    gname, best_match = gene_attrs_dict["gname"],gene_attrs_dict["best_match"]
    source, evidence = gene_attrs_dict['source'], gene_attrs_dict['evidence']
    pident, qcovs = gene_attrs_dict['pident'], gene_attrs_dict['qcovs']

    if term in ["transcript", "gene"]:
        uid = tuid
        best_match = trans_attrs["best_match"]
        source = trans_attrs["source"]
        evidence = trans_attrs["evidence"]
        pident = trans_attrs["pident"]
        qcovs = trans_attrs["qcovs"]
        tid = trans_attrs["uid"]

    if term == "gene":
        gid = trans_attrs["gid"]

    if term == "parent_change":
        gid = parent_gene
        tid = uid

    if term == "transcript_locus":
        source, evidence = "locus", "locus"
        gid = parent_gene
        tid = tuid

    modified = modify_transcript_attr(uid, tid, count, gid, gname, source, qcovs, pident, evidence,
                                      ann_term, seen_ids)

    collect[tuid] = extract_attr(tran_attr, "ID")
    seen_ids.add(collect[tuid])

    return modified


def annotate_existing(attr, count, trans_vals, term, parent_gene, seen_ids=set()):

    vals1 = copy.deepcopy(trans_vals)
    attrs_dict = get_attrs_dict(attr)
    collect = dict()

    for item in vals1:
        chrom, start, end, strand, current_attr, feat = parse_row(item)
        trans_attrs = get_attrs_dict(current_attr)
        tid = trans_attrs['uid']
        ann_gene = feat == "gene"
        ann_transcript = feat == "transcript"

        if ann_gene:
            gid = parse_attrs(current_attr, "ID")
            modified = modify_gene_attrs(gid, attrs_dict)

        elif ann_transcript:
            modified = annotate_transcript(parent_gene=parent_gene, term=term, trans_attr=current_attr,
                                           seen_ids=seen_ids, collect=collect, gene_attrs_dict=attrs_dict)

        else:
            pid = collect[tid]
            modified = modify_attr(current_attr, pid, count, attrs_dict['gname'], term)

        # Modify the attribute column
        item[8] = modified

    return vals1


def handle_gene_overlap(data, count, prev_gene="", seen=set()):
    # Get the transcript counts for this gene
    # Keeping the annotations from transcript itself, change the  parent to augustus gene

    # CASE 3A : transcript is annotated and gene is annotated -
    #  Nothing to do except correct transcript suffix numbers and assigning correct parent id
    if is_augustus(data=data):

        annotate = "transcript"
        target = data['gene_meta']['attr']
        parent_gene = parse_attrs(target, "ID")

        modified_transcript = annotate_existing(gene_attr, tcount, transcript_vals, annotate, parent_gene,
                                                seen_tids)
        if gene_key not in seen:
            merged_store.extend(aug_store[gene_key])
            seen.add(gene_key)
        merged_store.extend(modified_transcript)

        return


    return


def inc_count(current_key, prev_gene, count):
    inc = 1
    # This transcript has been seen before
    #TODO: this only works if it is immedtialy seen before.
    if current_key == prev_gene:
        inc = count + 1

    return inc


def create_merged_gff(gene_file, transcript_file):
    gstore, gmeta_data = parse_file(gene_file, "gene")
    tstore, tmeta_data = parse_file(transcript_file, "transcript")
    prev_gene, count = "", 0

    for tkey, tvals in tstore.items():
        # Get the transcript meta data and parse
        trans_meta = tmeta_data.get(tkey)
        data = parse_data(gmeta_data=gmeta_data, trans_meta=trans_meta, gstore=gstore)
        is_gene_overlap = data['is_gene_overlap']

        # CASE 1 : # transcript overlaps with more than one gene - remove transcript.
        if data['is_overlapped']:
            continue

        # CASE 2 : discard if the overlapped gene and transcript do not have the same orientation
        if is_gene_overlap and data['gene_strand'] != data['trans_strand']:
            continue

        # CASE 3 : Overlap with an augustus gene
        if is_gene_overlap:
            tcount = get_count(rows=data['gene'], target_str='transcript', target_idx=2)
            # Increment the transcript count depending
            # on how many times this transcript has already been seen
            count = inc_count(current_key=data['g_key'], prev_gene=prev_gene, count=count)
            tcount += count
            prev_gene = data['g_key']
            # Continue once the different types of gene overlaps are handled.
            if handle_gene_overlap(data=data, count=count, prev_gene=prev_gene):
                continue

        # CASE 4 : No overlap - add gene row to the largest transcript
        if data['no_overlap']:
            handle_no_overlap()
            continue



        1 / 0
        pass

    return


def main():

    create_merged_gff(gene_file=TEST_AUGUSTUS, transcript_file=TEST_STRINGTIE)

    return


if __name__ == "__main__":
    main()
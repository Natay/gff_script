import csv
import copy

TEST_AUGUSTUS = "aug.gff3"
TEST_STRINGTIE = "stringtie.gff3"


def parse_row(row):
    # chrom, start, end, strand, attr, feat
    return row[0], int(row[3]), int(row[4]), row[6], row[8], row[2]


def dict_to_str(data, delim=";"):
    """
    Convert dictionary to string or separated by delim.
    """
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


def modify_stringtie_gene_attr(gid, gene_id, gname, source, qcovs, pident,evidence):

    new_attr = dict(ID=gid, gene_id=gid, gene_name=gene_id,
                    Name=gene_id, best_match="stringtie")

    if not (evidence == "stringtie" or evidence == "augustus"):
        new_attr.update(dict(gene_name=gname, Name=gname, best_match=source,
                             query_coverage=qcovs, percentage_identity=pident))

    attr_str = dict_to_str(data=new_attr)
    return attr_str


def modify_transcript_attr(data):
    # TODO: working on this
    parent = gid
    trans_no = tid.split(".")[1]
    uid_mod = gid + ".t" + str(tidx)

    while uid_mod in seen_tuids:
        tidx += 1
        uid_mod = gid + ".t" + str(tidx)

    trans_id = gid + ".t" + str(tidx)
    attrs_dict = dict(ID=uid_mod, Parent=parent, original_id=uid,
                      transcript_id=trans_id,
                      gene_name=parent, Name=trans_id, best_match= evidence)

    if evidence == "stringtie":
        trans_id = gid + ".t" + str(tidx) if term != "gene" else uid_mod
        attrs_dict.update(dict(transcript_id=trans_id, Name=trans_id))

    elif evidence == "locus":
        trans_id = gname + ".t" + str(tidx) if term != "gene"  else ".".join([gname, trans_no])
        # gene_name={gname};Name={trans_id};best_match={evidence}
        attrs_dict.update(dict(transcript_id=trans_id, gene_name=gname, Name=trans_id, best_match=evidence))

    else:
        trans_id = gname + ".t" + str(tidx) if term != "gene"  else ".".join([gname, trans_no])

        if term == "stringtie_transcript" or term == "parent_change":
            attrs_dict.update(dict(transcript_id=trans_id, gene_name=gname,
                                   Name=trans_id,
                                   best_match=source,
                                   query_coverage=qcovs, percentage_identity=pident))
        else:
            attrs_dict.update(dict(ID=uid, transcript_id=trans_id, gene_name=gname,
                                   Name=trans_id,
                                   best_match=source,
                                   query_coverage=qcovs, percentage_identity=pident))
    new_attr = dict_to_str(data=attrs_dict)
    return new_attr


def modify_attr(attr, tid, tidx, gname, term):
    #TODO: working on this
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

    # TODO: ensure that seen_ids are actually changing correctly.
    data = dict(tid=tid, uid=uid, gid=gid, gname=gname, source=source,
                qcovs=qcovs, pident=pident, evidence=evidence, term=term, seen_ids=seen_ids)
    modified = modify_transcript_attr(data)
    collect[tuid] = parse_attrs(trans_attr, "ID")
    # TODO: ensure that seen_ids are actually changing correctly.
    seen_ids.add(collect[tuid])

    return modified


def annotate_existing(attr, count, trans_vals, annotate, parent_gene, seen_ids=set()):

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
            modified = annotate_transcript(parent_gene=parent_gene, term=annotate,
                                           trans_attr=current_attr,
                                           seen_ids=seen_ids,
                                           collect=collect, gene_attrs_dict=attrs_dict)
        else:
            pid = collect[tid]
            modified = modify_attr(current_attr, pid, count, attrs_dict['gname'], annotate)

        # Modify the attribute column
        item[8] = modified

    return vals1


def add_to_merge_store(merged_store, gene_vals, gene_key, trans_vals, seen=set()):

    if gene_key not in seen:
        merged_store.extend(gene_vals)
        seen.add(gene_key)
    merged_store.extend(trans_vals)


def find_gene_cords(tmeta, merged):
    """
    Find gene coordinates given a list of merged transcripts.
    """
    gstart, gend = tmeta['start'], tmeta['end']

    for item in merged:
        chrom, start, end, strand, tid = item
        start, end = int(start), int(end)
        # Construct new gene meta data with new values.
        local_gmeta = dict(chrom=chrom, start=start, end=end, strand=strand)

        if is_overlap(gmeta=local_gmeta, tmeta=tmeta):
            # Redefine gene start and end coordinates.
            gstart = tmeta['start'] if tmeta['start'] < start else start
            gend = tmeta['end'] if tmeta['end'] > end else end

    return gstart, gend


def handle_no_overlap(tmeta, tkey, tvals, merged_trans_list=[], gene_count=0, merged_gene=dict()):

    gchrom, gstart, gend, gstrand = tmeta['chrom'], tmeta['start'], tmeta['end'], tmeta['strand']
    tattrs = tmeta['attr']

    if merged_trans_list:
        gstart, gend = find_gene_cords(tmeta=tmeta,merged=merged_trans_list)

    best_match = parse_attrs(tattrs, "best_match")
    source = best_match
    tid = parse_attrs(tattrs, "ID")

    gid = ".".join(tid.split(".")[:-1])
    gname = parse_attrs(tattrs, "gene_name")

    evidence = get_evidence(best_match)
    if evidence == "augustus" or evidence == "stringtie":
        gname = gid
    pident, qcovs = get_evidence_vals(evidence, tattrs)
    #  There are transcripts which stringtie thinks should come from the same gene.
    # But augustus predicted two genes in that location instead of one which also got annotated.
    # So, assign uniq ids to these genes.
    stringtie_gene_id = gid
    gid = gid + ".gene" + str(gene_count)

    gene_attr = modify_stringtie_gene_attr(gid, stringtie_gene_id, gname, source, qcovs, pident, evidence)
    gene_row = [gchrom, "StringTie", "gene", gstart, gend, ".", gstrand, ".", gene_attr]

    if gname not in merged_gene:
        merged_gene[gname].append(gene_row)
    else:
        merged_gene[gname][0] = gene_row

    merged_gene.setdefault(gname, []).extend(tvals)
    merged_trans_list.append((gchrom, gstart, gend, gstrand, tkey))


def handle_gene_overlap(data, count, gstore, tvals, trans_attrs, merged_store=[], seen=set()):
    # TODO: working on this
    # Get the transcript counts for this gene
    # Keeping the annotations from transcript itself, change the  parent to augustus gene
    gene_key = data['g_key']
    gene_vals = gstore[gene_key]
    target = data['gene_meta'].get('attr')
    parent_gene = parse_attrs(target, "ID")

    not_aug = data['t_best'] != "stringtie" and data['g_best'] != "augustus"
    is_aug = data['t_best'] == "stringtie" and data['g_best'] == "augustus"
    cond3 = data['t_best'] != "stringtie" and data['g_best'] == "augustus"
    cond4 = data['t_best'] == "stringtie" and data['g_best'] != "augustus"

    # CASE 3A : transcript is annotated and gene is annotated -
    #  Nothing to do except correct transcript suffix numbers and assigning correct parent id
    if not_aug:
        annotate = "transcript"
        modified_transcript = annotate_existing(attr=target, count=count, trans_vals=tvals,
                                                annotate=annotate, parent_gene=parent_gene,
                                                seen_ids=seen)
        add_to_merge_store(merged_store=merged_store,
                           gene_vals=gene_vals,
                           gene_key=gene_key,
                           trans_vals=modified_transcript,
                           seen=seen)
        return True
    # CASE 3B : transcript is unannotated and gene is unannotated
    # Nothing to do except correct transcript suffix numbers and assign correct parent id
    if is_aug:
        annotate = "transcript_locus"
        modified_transcript = annotate_existing(attr=target, count=count, trans_vals=tvals,
                                                annotate=annotate, parent_gene=parent_gene,
                                                seen_ids=seen)
        add_to_merge_store(merged_store=merged_store,
                           gene_vals=gene_vals,
                           gene_key=gene_key,
                           trans_vals=modified_transcript,
                           seen=seen)
        return True

    # CASE 3C : transcript is annotated gene is unannotated - annotate gene with transcript ann
    # change transcript uid and parent
    if cond3:
        annotate = "gene"
        gene_attr = gene_vals[0][8]
        parent_gene = parse_attrs(gene_attr, "ID")
        modified_gene = annotate_existing(attr=trans_attrs, count=count, trans_vals=gene_vals,
                                          annotate=annotate, parent_gene=parent_gene,
                                          seen_ids=seen)

        modified_transcript = annotate_existing(attr=trans_attrs, count=count, trans_vals=tvals,
                                                annotate="parent_change", parent_gene=parent_gene,
                                                seen_ids=seen)
        add_to_merge_store(merged_store=merged_store,
                           gene_vals=modified_gene,
                           gene_key=gene_key,
                           trans_vals=modified_transcript,
                           seen=seen)
        return True

    # CASE 3D : transcript is unannotated and gene is annotated - annotate transcript with gene ann
    if cond4:
        annotate = "transcript_locus"
        gene_attr = gene_vals[0][8]
        parent_gene = parse_attrs(gene_attr, "ID")
        modified_transcript = annotate_existing(attr=gene_attr, count=count, trans_vals=tvals,
                                                annotate=annotate, parent_gene=parent_gene, seen_ids=seen)

        add_to_merge_store(merged_store=merged_store,
                           gene_vals=gene_vals,
                           gene_key=gene_key,
                           trans_vals=modified_transcript,
                           seen=seen)
        return True

    return False


def inc_count(current_key, prev_gene, count):
    inc = 1
    # This transcript has been seen before
    if current_key == prev_gene:
        inc = count + 1

    return inc


def create_merged_gff(gene_file, transcript_file):
    gstore, gmeta_data = parse_file(gene_file, "gene")
    tstore, tmeta_data = parse_file(transcript_file, "transcript")
    prev_gene, count = "", 0
    merged_trans, merged_gene = [], dict()
    gene_count = 0

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
            trans_attrs = trans_meta['attrs']
            # Continue once the different types of gene overlaps are handled.
            if handle_gene_overlap(data=data, count=tcount, tvals=tvals, gstore=gstore, trans_attrs=trans_attrs,
                                   merged_store=merged_trans):
                continue

        # CASE 4 : No overlap - add gene row to the largest transcript
        if data['no_overlap']:
            gene_count += 1
            if handle_no_overlap(merged_trans_list=merged_trans, tkey=tkey, tmeta=tmeta_data,
                                 tvals=tvals, merged_gene=merged_gene, gene_count=gene_count):
                continue

    # change the parent attribute of the strigtie-specific transcripts
    # to the newly assigned gene-id.
    # TODO: working on this
    new_gene_store = modify_transcript_parent(new_gene_store)
    new_gene_store = modify_stringtie_transcript_count(new_gene_store, seen_tids)
    return merged_store, new_gene_store



def process_output(merged_store, new_gene_store):
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

    return


def main():

    # Create merged gff file
    merged_store, new_gene_store = create_merged_gff(gene_file=TEST_AUGUSTUS, transcript_file=TEST_STRINGTIE)

    # Process the output.
    process_output(merged_store=merged_store, new_gene_store=new_gene_store)
    return


if __name__ == "__main__":
    main()
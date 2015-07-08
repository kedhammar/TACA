import re
import subprocess
import gzip
import glob
import os
import logging
import flowcell_parser.classes as cl
from taca.utils.config import CONFIG

logger=logging.getLogger(__name__)

dmux_folder='Demultiplexing'


def compute_undetermined_stats(run_dir, run_type, dex_status='COMPLETED'):
    """Will compute undetermined stats for all the undetermined files produced so far. 
    compute_index_freq takes care of concurrent runs.
    :param run_dir: path to the flowcell
    :type run_dir: str
    :param run_type: type of run (HiSeqX, HiSeq)
    :type run_type: string
    :param dex_status: status of the demux (COMPLETED/IN_PROGRESS)
    :type dex_status: string
    """
    global dmux_folder
    try:
        dmux_folder=CONFIG['analysis'][run_type]['bcl2fastq']['options'][0]['output_dir']
    except KeyError:
        dmux_folder='Demultiplexing'
    if os.path.exists(os.path.join(run_dir, dmux_folder)):
        xtp=cl.XTenParser(run_dir)
        ss=xtp.samplesheet
        workable_lanes=get_workable_lanes(run_dir, dex_status)
        for lane in workable_lanes:
            compute_index_freq(run_dir, lane)

    else:
        logger.warn("No demultiplexing folder found, aborting")



def check_lanes_QC(run, run_type,
                   max_percentage_undetermined_indexes_pooled_lane=5,
                   max_percentage_undetermined_indexes_unpooled_lane=20,
                   minimum_percentage_Q30_bases_per_lane=75,
                   minimum_yield_per_lane=305000000,
                   max_frequency_most_represented_und_index_pooled_lane=5,
                   max_frequency_most_represented_und_index_unpooled_lane=40,
                   dex_status='COMPLETED'):
    """Will check for undetermined fastq files, and perform the linking to the sample folder if the
    quality thresholds are met.

    :param run: path of the flowcell
    :type run: str
    :param max_percentage_undetermined_indexes_pooled_lane: max percentage of undetermined indexed allowed in a pooled lane
    :type max_percentage_undetermined_indexes_pooled_lane: float
    :param max_percentage_undetermined_indexes_unpooled_lane: max percentage of undetermined indexed allowed in a unpooled lane (one sample only)
    :type max_percentage_undetermined_indexes_unpooled_lane: float
    :param minimum_yield_per_lane: minimum lane yield
    :type minimum_yield_per_lane: int
    :param max_frequency_most_represented_und_index_pooled_lane: maximum frequency among the undetermined in a pooled lane
    :type max_frequency_most_represented_und_index_pooled_lane: float
    :param max_frequency_most_represented_und_index_unpooled_lane:  maximum frequency among the undetermined in a unpooled lane
    :type max_frequency_most_represented_und_index_unpooled_lane: float
    :param dex_status: status of the demultiplexing
    :type dex_status: sting (COMPLETED/RUNNING)

    :returns boolean: True  if the flowcell passes the checks, False otherwise
    """
    

    status=True
    if os.path.exists(os.path.join(run, dmux_folder)):
        xtp=cl.XTenParser(run)
        
        global dmux_folder
        try:
            dmux_folder=CONFIG['analysis'][run_type]['bcl2fastq']['options'][0]['output-dir']
        except KeyError:
            dmux_folder='Demultiplexing'

        
        ss=xtp.samplesheet
        lb=xtp.lanebarcodes
        lanes=xtp.lanes
        if not lb:
            logger.info("The HTML report is not available. QC cannot be performed, the FC will not be tranferred.")
            return False


        #these two functions need to became lane specific and return an array
        path_per_lane=get_path_per_lane(run, ss)
        samples_per_lane=get_samples_per_lane(ss)
        workable_lanes=get_workable_lanes(run, dex_status)
        for lane in workable_lanes:
            #QC on the yield
            lane_status = True
            if lane_check_yield(lane, lanes, minimum_yield_per_lane):
                lane_status = lane_status and True
            else:
                logger.warn("lane {} did not pass yield qc check. This FC will not be transferred.".format(lane))
                lane_status = lane_status and False
            #QC on the total %>Q30 of the all lane
            if lane_check_Q30(lane, lanes, minimum_percentage_Q30_bases_per_lane):
                lane_status = lane_status and True
            else:
                logger.warn("lane {} did not pass Q30 qc check. This FC will not be transferred.".format(lane))
                lane_status = lane_status and False
            
            #now QC for undetermined
            max_percentage_undetermined_indexes = max_percentage_undetermined_indexes_pooled_lane
            max_frequency_most_represented_und  = max_frequency_most_represented_und_index_pooled_lane
            #distinguish the case between Pooled and Unpooled lanes, for unpooled lanes rename the Undetemriend file
            if is_unpooled_lane(ss, lane):
                rename_undet(run, lane, samples_per_lane)
                max_percentage_undetermined_indexes = max_percentage_undetermined_indexes_unpooled_lane
                max_frequency_most_represented_und  = max_frequency_most_represented_und_index_unpooled_lane

            if check_undetermined_reads(run, lane, lanes, max_percentage_undetermined_indexes):
                if check_maximum_undertemined_freq(run, lane, max_frequency_most_represented_und):
                    if is_unpooled_lane(ss, lane):
                        link_undet_to_sample(run, lane, path_per_lane)
                    lane_status= lane_status and True
                else:
                    logger.warn("lane {} did not pass the undetermiend qc checks. Frequency of most popular undetermined index is too large.".format(lane))
                    lane_status= lane_status and False
            else:
                logger.warn("lane {} did not pass the undetermiend qc checks. Fraction of undetermined too large.".format(lane))
                lane_status= lane_status and False
            if lane_status:
                logger.info("lane {} passed all qc checks".format(lane))
            #store the status for the all FC
            status = status and lane_status

    else:
        logger.warn("No demultiplexing folder found, aborting")
        status = False
    
    return status
        
        

def qc_for_pooled_lane(lane,lb , und_thresh):
    """Performs QC checks for pooled lanes, i.e., more than one sample per lane
        
    :param lane: lane number (form 1 to 8)
    :type lane: int
    :param lb: parsed description of the lane
    :type lb: XTenParse.lanebarcodes object
    :param und_thresh: maximum allowed percentage of undetemined reads
    :type und_thresh: int
    """
    d={}
    d['det']=0
    d['undet']=0
    for entry in lb.sample_data:
        if lane == int(entry['Lane']):
            if entry.get('Sample')!='unknown':
                d['det']+=int(entry['Clusters'].replace(',',''))
            else:
                d['undet']=int(entry['Clusters'].replace(',',''))


    if d['undet'] > (d['det']+d['undet']) * und_thresh / 100:
        logger.warn("Lane {} has more than {}% undetermined indexes ({}%)".format(lane, und_thresh,d['undet']/(d['det']+d['undet'])*100))
        return False

    return True

    


def rename_undet(run, lane, samples_per_lane):
    """Renames the Undetermined fastq file by prepending the sample name in front of it

    :param run: the path to the run folder
    :type run: str
    :param status: the demultiplexing status
    :type status: str
    :param samples_per_lane: lane:sample dict
    :type status: dict
    """
    for file in glob.glob(os.path.join(run, dmux_folder, "Undetermined*L0?{}*".format(lane))):
        old_name=os.path.basename(file)
        old_name_comps=old_name.split("_")
        old_name_comps[1]=old_name_comps[0]# replace S0 with Undetermined
        old_name_comps[0]=samples_per_lane[lane]#replace Undetermined with samplename
        for index, comp in enumerate(old_name_comps):
            if comp.startswith('L00'):
                old_name_comps[index]=comp.replace('L00','L01')#adds a 1 as the second lane number in order to differentiate undetermined from normal in piper
                
        new_name="_".join(old_name_comps)
        logger.info("Renaming {} to {}".format(file, os.path.join(os.path.dirname(file), new_name)))
        os.rename(file, os.path.join(os.path.dirname(file), new_name))


def get_workable_lanes(run, status):
    """List the lanes that have a .fastq file

    :param run: the path to the run folder
    :type run: str
    :param status: the demultiplexing status
    :type status: str

    :rtype: list of ints 
    :returns:: list of lanes having an undetermined fastq file
    """
    lanes=[]
    pattern=re.compile('L0[0,1]([0-9])')
    for unde in glob.glob(os.path.join(run, dmux_folder, '*Undetermined_*')):
        name=os.path.basename(unde)
        lanes.append(int(pattern.search(name).group(1)))
    lanes=list(set(lanes))
    if status =='IN_PROGRESS': 
        #the last lane is the one that is currently being worked on by bcl2fastq, don't work on it.
        lanes=lanes[:-1]
    logger.info("post_demux processing will happen with lanes {}".format(lanes))
    return lanes


def link_undet_to_sample(run, lane, path_per_lane):
    """symlinks the undetermined file to the right sample folder with a RELATIVE path so it's carried over by rsync
    
    :param run: path of the flowcell
    :type run: str
    :param lane: lane identifier
    :type lane: int
    :param path_per_lane: {lane:path/to/the/sample}
    :type path_per_lane: dict"""
    for fastqfile in glob.glob(os.path.join(run, dmux_folder, '*Undetermined*_L0?{}_*'.format(lane))):
        if not os.path.exists(os.path.join(path_per_lane[lane], os.path.basename(fastqfile))):
            fqbname=os.path.basename(fastqfile)
            logger.info("linking file {} to {}".format(fastqfile, path_per_lane[lane]))
            os.symlink(os.path.join('..','..',fqbname), os.path.join(path_per_lane[lane], os.path.basename(fastqfile)))

def save_index_count(barcodes, run, lane):
    """writes the barcode counts

    :param barcodes: {barcode:count}
    :type barcodes: dict
    :param run: path to the flowcell
    :type run: str
    """
    with open(os.path.join(run, dmux_folder, 'index_count_L{}.tsv'.format(lane)), 'w') as f:
        for barcode in sorted(barcodes, key=barcodes.get, reverse=True):
            f.write("{}\t{}\n".format(barcode, barcodes[barcode]))


def compute_index_freq(run, lane):
    """uses subprocess to perform zcat <file> | sed -n '1~4 p' | awk -F ':' '{print $NF}', counts the barcodes and 
    returns true if the most represented index accounts for less than freq_tresh% of the total
    The functions creates a semaphor file, .running that prevents cuncurrent starts. This flag file(s) is used to avoid
    early starts of the final step due to concurrent starts.
    
    :param run: path to the flowcell
    :type run: str
    :param lane: lane identifier
    :type lane: int
    """
    barcodes={}
    if os.path.exists(os.path.join(run, dmux_folder,'index_count_L{}.running'.format(lane))):
        logger.info("Lane {} is currenlty under processing: WARNING this lane is taking too long.".format(lane))
        return
    elif os.path.exists(os.path.join(run, dmux_folder,'index_count_L{}.tsv'.format(lane))):
        logger.info("Found index count for lane {}. No need to recompute it.".format(lane))
        return
    else:
        open(os.path.join(run, dmux_folder,'index_count_L{}.running'.format(lane)), 'a').close()
        logger.info("Lane {} is currenlty under processing: counting undetermined indexes.".format(lane))
    #next check is not needed but I want to be verbose
    if os.path.exists(os.path.join(run, dmux_folder,'index_count_L{}.tsv'.format(lane))):
        logger.info("Lane {} raised an exception: .tsv file present but should not.".format(lane))
        return
    open(os.path.join(run, dmux_folder,'index_count_L{}.tsv'.format(lane)), 'a').close()
    for fastqfile in glob.glob(os.path.join(run, dmux_folder, 'Undetermined_S*_L00{}_R1*'.format(lane))):
            logger.info("working on {}".format(fastqfile))
            zcat=subprocess.Popen(['zcat', fastqfile], stdout=subprocess.PIPE)
            sed=subprocess.Popen(['sed', '-n', "1~4p"],stdout=subprocess.PIPE, stdin=zcat.stdout)
            awk=subprocess.Popen(['awk', '-F', ":", '{print $NF}'],stdout=subprocess.PIPE, stdin=sed.stdout)
            zcat.stdout.close()
            sed.stdout.close()
            output = awk.communicate()[0]
            zcat.wait()
            sed.wait()
            for barcode in output.split('\n')[:-1]:
                try:
                    barcodes[barcode]=barcodes[barcode]+1
                except KeyError:
                    barcodes[barcode]=1
    save_index_count(barcodes, run, lane)
    #now remove the flag file
    os.remove(os.path.join(run, dmux_folder,'index_count_L{}.running'.format(lane)))
    return


def check_maximum_undertemined_freq(run, lane, freq_tresh):
    """uses the counts made by compute_index_freq to  counts the barcodes and
        returns true if the most represented index accounts for less than freq_tresh% 
        of the total amount of undetermiend
        
        :param run: path to the flowcell
        :type run: str
        :param lane: lane identifier
        :type lane: int
        :param freq_tresh: maximal allowed frequency of the most frequent undetermined index
        :type frew_tresh: float
        :rtype: boolean
        :returns: True if the checks passes, False otherwise
        """
    barcodes={}
    if os.path.exists(os.path.join(run, dmux_folder,'index_count_L{}.tsv'.format(lane))):
        with open(os.path.join(run, dmux_folder,'index_count_L{}.tsv'.format(lane))) as idxf:
            for line in idxf:
                barcodes[line.split('\t')[0]]=int(line.split('\t')[1])
    else:
        logger.warn("file index_count_L{}.tsv not found. This FC will be marked as QC failed and not tranfrred".format(lane))
        return False

    total=sum(barcodes.values())
    count, bar = max((v, k) for k, v in barcodes.items())
    if total * (freq_tresh / float(100)) < count:
        logger.warn("The most frequent barcode of lane {} ({}) represents {}%, "
                    "which is over the threshold of {}%".format(lane, bar, count * 100 / total , freq_tresh))
        return False
    else:
        return True

                    

def check_undetermined_reads(run, lane, lanes, freq_tresh):
    """uses the counts made by compute_index_freq to  counts the barcodes and
    returns true if the most represented index accounts for less than freq_tresh% of the total
                
    :param run: path to the flowcell
    :type run: str
    :param lane: lane identifier
    :type lane: int
    :param lanes: reader of lanes.html
    :type lanes: flowcell_parser.classes.XTenLaneBarcodeParser
    :param freq_tresh: maximal allowed percentage of undetermined indexes in a lane
    :type frew_tresh: float
    :rtype: boolean
    :returns: True if the checks passes, False otherwise
    """
    #compute lane yield
    lane_yield = 0;
    for entry in lanes.sample_data:
        if lane == int(entry['Lane']):
            if lane_yield > 0:
                logger.warn("lane_yeld must be 0, somehting wrong is going on here")
            lane_yield = int(entry['Clusters'].replace(',',''))

    barcodes={}
    if os.path.exists(os.path.join(run, dmux_folder,'index_count_L{}.tsv'.format(lane))):
        with open(os.path.join(run, dmux_folder,'index_count_L{}.tsv'.format(lane))) as idxf:
            for line in idxf:
                    barcodes[line.split('\t')[0]]=int(line.split('\t')[1])
    else:
        logger.warn("file index_count_L{}.tsv not found. This FC will be marked as QC failed and not tranfrred".format(lane))
        return False


    undetermined_total=sum(barcodes.values())
    percentage_und = (undetermined_total/float(lane_yield))*100
    if  percentage_und > freq_tresh:
        logger.warn("The undetermined indexes account for {}% of lane {}, "
                    "which is over the threshold of {}%".format(percentage_und, lane, freq_tresh))
        return False
    else:
        return True
                    


def lane_check_yield(lane, lanes, minimum_yield):
    """checks that the total yield lane (P/F reads) is higher than the minimum
    :param lane: lane currenlty being worked
    :type lane: int
    :param lanes: reader of lanes.html
    :type lanes: flowcell_parser.classes.XTenLaneBarcodeParser
    :param minimum_yield: minimum yield as specified by documentation
    :type minimum_yield: float
    
    :rtype: boolean
    :returns: True if the lane has an yield above the specified minimum
    """
    
    for entry in lanes.sample_data:
        if lane == int(entry['Lane']):
            lane_clusters = int(entry['Clusters'].replace(',',''))
            if lane_clusters >= minimum_yield:
                return True
    return False


def lane_check_Q30(lane, lanes, q30_tresh):
    """checks that the total Q30 of the lane  is higher than the minimum
    :param lane: lane currenlty being worked
    :type lane: int
    :param lanes: reader of lanes.html
    :type lanes: flowcell_parser.classes.XTenLaneBarcodeParser
    :param q30_tresh: Q30 threshold
    :type q30_tresh: float
        
    :rtype: boolean
    :returns: True if the lane has a Q30 above the specified minimum
    """

    for entry in lanes.sample_data:
        if lane == int(entry['Lane']):
            if float(entry['% >= Q30bases']) >= q30_tresh:
                return True
    return False






def sample_qc_check(lane, lb, und_tresh, q30_tresh):
    """checks wether the percentage of bases over q30 for the sample is
    above the treshold, and if the amount of undetermined is below the treshold
    
    :param lane: lane identifier
    :type lane: int
    :param lb: reader of laneBarcodes.html
    :type lb: flowcell_parser.classes.XTenLaneBarcodes
    :param und_tresh: maximal allowed percentage of undetermined indexes
    :type und_tresh: float
    :param q30_tresh: maximal allowed  percentage of bases over q30
    :type q30_tresh: float

    :rtype: boolean
    :returns: True of the qc checks pass, False otherwise
    
    """
    d={}
    for entry in lb.sample_data:
        if lane == int(entry['Lane']):
            if entry.get('Sample')!='unknown':
                if float(entry['% >= Q30bases']) < q30_tresh:
                    logger.warn("Sample {} of lane {} has a percentage of bases over q30 of {}%, "
                            "which is below the cutoff of {}% ".format(entry['Sample'], lane, float(entry['% >= Q30bases']), q30_tresh))
                    return False
                d['det']=int(entry['Clusters'].replace(',',''))
            else:
                d['undet']=int(entry['Clusters'].replace(',',''))


    if d['undet'] > (d['det']+d['undet']) * und_tresh / 100:
        logger.warn("Lane {} has more than {}% undetermined indexes ({:.2f}%)".format(lane, und_tresh,float(d['undet'])/(d['det']+d['undet'])*100))
        return False

    return True




def get_path_per_lane(run, ss):
    """
    :param run: the path to the flowcell
    :type run: str
    :param ss: SampleSheet reader
    :type ss: flowcell_parser.XTenSampleSheet
    """
    d={}
    for l in ss.data:
        try:
            d[int(l['Lane'])]=os.path.join(run, dmux_folder, l['Project'], l['SampleID'])
        except KeyError:
            logger.error("Can't find the path to the sample, is 'Project' in the samplesheet ?")
            d[int(l['Lane'])]=os.path.join(run, dmux_folder)

    return d

def get_samples_per_lane(ss):
    """
    :param ss: SampleSheet reader
    :type ss: flowcell_parser.XTenSampleSheet
    :rtype: dict
    :returns: dictionnary of lane:samplename
    """
    d={}
    for l in ss.data:
        s=l['SampleName'].replace("Sample_", "").replace("-", "_")
        d[int(l['Lane'])]=l['SampleName']

    return d

def get_barcode_per_lane(ss):
    """
    :param ss: SampleSheet reader
    :type ss: flowcell_parser.XTenSampleSheet
    :rtype: dict
    :returns: dictionnary of lane:barcode
    """
    d={}
    for l in ss.data:
        d[int(l['Lane'])]=l['index']

    return d


def is_unpooled_lane(ss, lane):
    """
    :param ss: SampleSheet reader
    :type ss: flowcell_parser.XTenSampleSheet
    :param lane: lane identifier
    :type lane: int
    :rtype: boolean
    :returns: True if the samplesheet has one entry for that lane, False otherwise
    """
    count=0
    for l in ss.data:
        if int(l['Lane']) == lane:
            count+=1
    return count==1

def is_unpooled_run(ss):
    """
    :param ss: SampleSheet reader
    :type ss: flowcell_parser.XTenSampleSheet
    :rtype: boolean
    :returns: True if the samplesheet has one entry per lane, False otherwise
    """
    ar=[]
    for l in ss.data:
        ar.append(l['Lane'])
    return len(ar)==len(set(ar))

        

if __name__=="__main__":
    import sys

    mainlog = logging.getLogger(__name__)
    mainlog.setLevel(level=logging.INFO)
    mfh = logging.StreamHandler(sys.stderr)
    mft = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    mfh.setFormatter(mft)
    mainlog.addHandler(mfh)
        
    check_undetermined_status("/srv/illumina/HiSeq_X_data/150424_ST-E00214_0031_BH2WY7CCXX")

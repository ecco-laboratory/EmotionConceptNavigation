#----------------------ASK THIS!---------------------------------------------------
#Qs: Does run1 and run2 naming depend on stimuli or order? 
# #Do i do odd subs, face_short run1, face_long run2 and even subs, word_short run1, word_long run2? 
#---------------------------------------------------------------------------

#template heuristic
#notes by EC Hahn, some additions by Monica Thieu

# Internal function, you shouldn't need to touch this
def create_key(template, outtype=('nii.gz',), annotation_classes=None):
    if template is None or not template:
        raise ValueError('Template must be a valid format string')
    return template, outtype, annotation_classes

# If there is ambiguity about which fieldmap scans correspond to which functional scans,
# these options will tell heudiconv (and later fmriprep) how to assign them
# Delete 'CustomAcquisitionLabel' from matching_parameters if you are not using
# custom acquisition BIDS tags in your scan names like "acq-mb8", "acq-me", etc.
# Refer to https://heudiconv.readthedocs.io/en/latest/heuristics.html for other options
# Delete this if you don't have any fieldmap scans
POPULATE_INTENDED_FOR_OPTS = {
    'matching_parameters': ['CustomAcquisitionLabel', 'ImagingVolume', 'Shims'],
    'criterion': 'Closest'
}

def infotodict(seqinfo):
    """Heuristic evaluator for determining which runs belong where
    allowed template fields - follow python string module:
    item: index within category
    subject: participant id
    seqitem: run number during scanning
    subindex: sub index within group
    """

    stem_subject = 'sub-{subject}'
    if_subject_odd = False if int(stem_subject.split('-')[1]) % 2 == 0 else True
    faceshort_run = 'run-01' if if_subject_odd else 'run-02'
    facelong_run = 'run-02' if if_subject_odd else 'run-01'

    #here is where you start changing protocol names to match what is in the dicom summary file. 
    #Each task needs to have formatting that keeps the first "task". 
    #For example, for the first run of the bandit task, the last part of the line 
    #should read "task-Bandit1-1_bold"

    stem_face = 'task-face'
    stem_word = 'task-word'
    stem_story = 'task-story'

    stem_anat = stem_subject + '/anat/'
    stem_fmap = stem_subject + '/fmap/'
    stem_func = stem_subject + '/func/'

    # the names in this dict should be partial unique string-matches
    # for the associated protocol names
    # it is important that each paired SBRef key be defined before its main sequence key
    # because of the string matching that will be used below
    run_types = {
        #keep these three the same
        't1_mpr': create_key(stem_anat + '_'.join([stem_subject, 'T1w'])), 
        'AP': create_key(stem_fmap + '_'.join([stem_subject, 'dir-AP', 'epi'])), 
        'PA': create_key(stem_fmap + '_'.join([stem_subject, 'dir-PA', 'epi'])),  
        'faceshort': create_key(stem_func + '_'.join([stem_subject, stem_face, faceshort_run, 'bold'])), 
        'facelong': create_key(stem_func + '_'.join([stem_subject, stem_face, facelong_run, 'bold'])),
        'word1': create_key(stem_func + '_'.join([stem_subject, stem_word, 'run-01', 'bold'])), 
        'word2': create_key(stem_func + '_'.join([stem_subject, stem_word, 'run-02', 'bold'])) 
        #'story': create_key(stem_func + '_'.join([stem_subject, stem_story, 'bold']))
    }

    # create a dict `info` such that each key has a separate empty list associated with it
    # tuples can be dictionary keys, believe it or not
    # but this double-dict system allows us to keep `info` in the format heudiconv expects
    # while keying in with a dict whose key names are more human-readable
    info = {key: [] for key in run_types.values()}
 
#change the dimensions (dim1, dim2, etc) and protocol names below to match the info in the dicom file
    for _, s in enumerate(seqinfo):

        for run_type in run_types.keys():
            if run_type in s.protocol_name:

                info[run_types[run_type]] = [s.series_id]
                break
     
    return info


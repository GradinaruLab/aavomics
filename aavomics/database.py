import os
import pandas
import numpy
import datetime

from enum import Enum

DATA_PATH = "data"

CELL_SETS = []
CELL_SETS_DICT = {}
TISSUE_SAMPLES = []
TISSUE_SAMPLES_DICT = {}
READ_SETS = []
READ_SETS_DICT = {}
SEQUENCING_LIBRARIES = []
SEQUENCING_LIBRARIES_DICT = {}
REFERENCES = []
REFERENCES_DICT = {}
ALIGNMENTS = []
ALIGNMENTS_DICT = {}
ANIMALS = []
ANIMALS_DICT = {}
TEMPLATES = []
TEMPLATES_DICT = {}
VIRUSES = []
VIRUSES_DICT = {}
VECTORS = []
VECTORS_DICT = {}
VECTOR_POOLS = []
VECTOR_POOLS_DICT = {}
INJECTIONS = []
INJECTIONS_DICT = {}
DISSOCIATION_RUNS = []
DISSOCIATION_RUNS_DICT = {}
SEQUENCING_RUNS = []
SEQUENCING_RUNS_DICT = {}

if "AAVOMICS_PATH" in os.environ:
    DATA_PATH = os.environ["AAVOMICS_PATH"]
    

class Reference():
    
    def __init__(self):
        self.name = ""

class Alignment():
    
    def __init__(self):
        self.name = ""
        self.reference = None
        self.aligner = ""
        self.aligner_version = 0

class Tissue_Sample():
    
    def __init__(self):
        self.name = ""
        self.region = None
        self.animal = None
        self.dissociation_run = None

class Cell_Set():

    def __init__(self):
        self.name = ""
        self.source_tissue = ""
        self.target_num_cells = 0
        self.sequencing_libraries = []
        
    def has_reads(self, read_type="Transcriptome"):
        
        for sequencing_library in self.sequencing_libraries:
            if sequencing_library.type != read_type:
                continue
            if len(sequencing_library.read_sets) > 0:
                return True
        
        return False
        
    def get_transcripts_directory(self, alignment_name=None, transcript_type="transcriptome"):
        
        transcripts_path = os.path.join(DATA_PATH, "cell_sets", self.name, transcript_type, "transcripts")
        
        # If an alignment is not specified, grab the first one
        if alignment_name is None:
            alignment_name = os.listdir(transcripts_path)[0]
            if len(os.listdir(transcripts_path)) > 1:
                print("WARNING: Multiple alignments found, picking %s." % alignment_name)
        
        transcripts_path = os.path.join(transcripts_path, alignment_name)
        
        if not os.path.exists(transcripts_path):
            os.makedirs(transcripts_path)
        
        return transcripts_path
        
    def get_alignment_directory(self, alignment_name=None, transcript_type="virus"):
        
        alignment_path = os.path.join(DATA_PATH, "cell_sets", self.name, transcript_type, "alignments")
        
        # If an alignment is not specified, grab the first one
        if alignment_name is None:
            alignment_name = os.listdir(alignment_path)[0]
            if len(os.listdir(alignment_path)) > 1:
                print("WARNING: Multiple alignments found, picking %s." % alignment_name)
        
        alignment_path = os.path.join(alignment_path, alignment_name)
        
        if not os.path.exists(alignment_path):
            os.makedirs(alignment_path)
        
        return alignment_path
        
    def get_virus_transcripts_directory(self, alignment_name=None):
        
        print("Deprecated, use get_transcripts_directory")
        return self.get_transcripts_directory(alignment_name=alignment_name, transcript_type="virus")
    
    def get_cell_ranger_h5_file_path(self, filtered=False, alignment_name=None):
        
        transcripts_path = self.get_transcripts_directory(alignment_name=alignment_name)
        
        if filtered:
            return os.path.join(transcripts_path, "filtered_feature_bc_matrix.h5")
        else:
            return os.path.join(transcripts_path, "raw_feature_bc_matrix.h5")
        
    def get_anndata_file_path(self, prefix="raw", alignment_name=None, transcript_type="transcriptome"):
        
        transcripts_path = self.get_transcripts_directory(alignment_name=alignment_name, transcript_type=transcript_type)
        
        return os.path.join(transcripts_path, "%s_barcode_transcript_counts.h5ad" % prefix)
        
class Sequencing_Library():
    
    def __init__(self):
        self.name = ""
        self.type = ""
        self.cell_set = None
        self.read_sets = []
        
class Read_Set():
    
    def __init__(self):
        self.name = ""
        self.sequencing_libraries = []
        self.sequencing_run = None
        
class Sequencing_Run():
    
    def __init__(self):
        self.name = ""
        self.read_1_length = 0
        self.read_2_length = 0
        self.read_sets = []

class Animal():
    
    def __init__(self):
        self.name = ""
        self.injections = []
        self.extraction_date = None
        self.DOB = None
        self.injection_date = None

class Injection():
    
    def __init__(self):
        self.animal = None
        self.vector_pool = None
        self.quantity = 0
        
class Vector_Pool():
    
    def __init__(self):
        self.id = ""
        self.name = ""
        self.vectors = []
        self.injections = []

class Vector():
    
    def __init__(self):
        self.id = ""
        self.name = ""
        self.delivery_vehicle = None
        self.cargo = None
        self.vector_pools = []
        self.barcode = None
        self.barcode_location = 0

class Template():
    
    def __init__(self):
        self.name = ""
        self.sequence = ""
        self.vectors = []

class Virus():
    
    def __init__(self):
        self.id = ""
        self.name = ""
        self.vectors = []

class Dissociation_Run():
    
    def __init__(self):
        self.name = ""
        self.tissue_samples = []
        self.protocol_version = 0

    
def load_database(data_path=None):
        
    global CELL_SETS
    global CELL_SETS_DICT
    global TISSUE_SAMPLES
    global TISSUE_SAMPLES_DICT
    global SEQUENCING_LIBRARIES
    global SEQUENCING_LIBRARIES_DICT
    global READ_SETS
    global READ_SETS_DICT
    global REFERENCES
    global REFERENCES_DICT
    global ALIGNMENTS
    global ALIGNMENTS_DICT
    global TEMPLATES
    global TEMPLATES_DICT
    global VIRUSES
    global VIRUSES_DICT
    global VECTORS
    global VECTORS_DICT
    global ANIMALS
    global ANIMALS_DICT
    global VECTOR_POOLS
    global VECTOR_POOLS_DICT
    global INJECTIONS
    global INJECTIONS_DICT
    global DISSOCIATION_RUNS
    global DISSOCIATION_RUNS_DICT
    global SEQUENCING_RUNS
    global SEQUENCING_RUNS_DICT
    
    if not data_path:
        data_path = DATA_PATH
        
    database_path = os.path.join(data_path, "database")
    database_file_paths = {}
    
    for file_name in os.listdir(database_path):
        data_type = ".".join(file_name.split(".")[0:-1])
        if data_type not in database_file_paths:
            database_file_paths[data_type] = os.path.join(database_path, file_name)
        else:
            print("Multiple database files found for %s. Using %s; ignoring others." % (data_type, database_file_paths[data_type]))
            
    tissue_samples_df = read_tissue_samples_df(database_file_paths["Tissue Samples"])
    
    for tissue_sample_row in tissue_samples_df.iterrows():
        
        tissue_sample = Tissue_Sample()
        tissue_sample.name = tissue_sample_row[0]
        tissue_sample.region = tissue_sample_row[1]["Region"]
        
        TISSUE_SAMPLES.append(tissue_sample)
        TISSUE_SAMPLES_DICT[tissue_sample.name] = tissue_sample
    
    cell_sets_df = read_cell_sets_df(database_file_paths["Cell Sets"])
    
    for cell_set_row in cell_sets_df.iterrows():
        
        cell_set = Cell_Set()
        cell_set.name = cell_set_row[0]
        cell_set.source_tissue = TISSUE_SAMPLES_DICT[cell_set_row[1]["Source Tissue"]]
        
        if not numpy.isnan(cell_set_row[1]["Target # Cells"]):
            cell_set.target_num_cells = int(cell_set_row[1]["Target # Cells"])
        
        CELL_SETS.append(cell_set)
        CELL_SETS_DICT[cell_set.name] = cell_set
    
    sequencing_libraries_df = read_sequencing_libraries_df(database_file_paths["Sequencing Libraries"])
    
    for sequencing_library_row in sequencing_libraries_df.iterrows():
        
        sequencing_library = Sequencing_Library()
        sequencing_library.name = sequencing_library_row[0]
        sequencing_library.type = sequencing_library_row[1]["Read 2 Type"]
        
        cell_set_name = sequencing_library_row[1]["Cell Set"]
        cell_set = CELL_SETS_DICT[cell_set_name]
        sequencing_library.cell_set = cell_set
        cell_set.sequencing_libraries.append(sequencing_library)
        
        SEQUENCING_LIBRARIES.append(sequencing_library)
        SEQUENCING_LIBRARIES_DICT[sequencing_library.name] = sequencing_library
        
    read_sets_df = read_read_sets_df(database_file_paths["Read Sets"])
    
    for read_set_row in read_sets_df.iterrows():
        
        read_set = Read_Set()
        read_set.name = read_set_row[0]
        read_set.num_reads = int(read_set_row[1]["# Reads"])
        
        for sequencing_library_name in read_set_row[1]["Libraries"].split(","):
            sequencing_library = SEQUENCING_LIBRARIES_DICT[sequencing_library_name]
            read_set.sequencing_libraries.append(sequencing_library)
            sequencing_library.read_sets.append(read_set)
            
        READ_SETS.append(read_set)
        READ_SETS_DICT[read_set.name] = read_set
        
    references_df = read_references_df(database_file_paths["References"])
    
    for reference_row in references_df.iterrows():
        
        reference = Reference()
        reference.name = reference_row[0]
        
        REFERENCES.append(reference)
        REFERENCES_DICT[reference.name] = reference
        
    alignments_df = read_alignments_df(database_file_paths["Alignments"])
    
    for alignment_row in alignments_df.iterrows():
        
        alignment = Alignment()
        alignment.name = alignment_row[0]
        alignment.aligner = alignment_row[1]["Aligner"]
        alignment.aligner_version = alignment_row[1]["Aligner Version"]
        
        reference_name = alignment_row[1]["Reference"]
        alignment.reference = REFERENCES_DICT[reference_name]
            
        ALIGNMENTS.append(alignment)
        ALIGNMENTS_DICT[alignment.name] = alignment
        
    templates_df = read_templates_df(database_file_paths["Templates"])
    
    for template_row in templates_df.iterrows():
        
        template = Template()
        template.name = template_row[0]
        template.sequence = template_row[1]["Sequence"]
        
        TEMPLATES.append(template)
        TEMPLATES_DICT[template.name] = template
        
    viruses_df = read_viruses_df(database_file_paths["Viruses"])
    
    for virus_row in viruses_df.iterrows():
        
        virus = Virus()
        virus.id = virus_row[0]
        virus.name = virus_row[1]["Name"]
        
        VIRUSES.append(virus)
        VIRUSES_DICT[virus.id] = virus
        
    vectors_df = read_vectors_df(database_file_paths["Vectors"])
    
    for vector_row in vectors_df.iterrows():
        
        vector = Vector()
        vector.id = vector_row[0]
        vector.name = vector_row[1]["Name"]
        virus = VIRUSES_DICT[vector_row[1]["Delivery Vehicle"]]
        vector.delivery_vehicle = virus
        virus.vectors.append(vector)
        cargo = TEMPLATES_DICT[vector_row[1]["Cargo"]]
        cargo.vectors.append(vector)
        vector.cargo = cargo
        if vector_row[1]["Barcode"]:
            vector.barcode = vector_row[1]["Barcode"]
        if vector_row[1]["Barcode Insertion Location"]:
            vector.barcode_location = vector_row[1]["Barcode Insertion Location"]
        
        VECTORS.append(vector)
        VECTORS_DICT[vector.id] = vector
        
    animals_df = read_animals_df(database_file_paths["Animals"])
    
    for animal_row in animals_df.iterrows():
        
        animal = Animal()
        animal.name = animal_row[0]
        
        if animal_row[1]["Samples"] is not None:
            for tissue_sample_name in animal_row[1]["Samples"].split(","):
                tissue_sample = TISSUE_SAMPLES_DICT[tissue_sample_name]
                tissue_sample.animal = animal
        if animal_row[1]["Tissue extraction date"] is not None and not isinstance(animal_row[1]["Tissue extraction date"], float):
            animal.extraction_date = datetime.datetime.strptime(animal_row[1]["Tissue extraction date"], "%b %d, %Y")
        if animal_row[1]["DOB"] is not None and not isinstance(animal_row[1]["DOB"], float):
            animal.DOB = datetime.datetime.strptime(animal_row[1]["DOB"], "%b %d, %Y")
        if animal_row[1]["Injection date"] is not None and not isinstance(animal_row[1]["Injection date"], float):
            animal.injection_date = datetime.datetime.strptime(animal_row[1]["Injection date"], "%b %d, %Y")
        
        ANIMALS.append(animal)
        ANIMALS_DICT[animal.name] = animal
        
    vector_pools_df = read_vector_pools_df(database_file_paths["Vector Pools"])
    
    for vector_pool_row in vector_pools_df.iterrows():
        
        vector_pool = Vector_Pool()
        vector_pool.id = vector_pool_row[0]
        vector_pool.name = vector_pool_row[1]["Name"]
        
        if vector_pool_row[1]["Vectors"] is not None:
            for vector_id in vector_pool_row[1]["Vectors"].split(","):

                vector = VECTORS_DICT[vector_id]
                vector_pool.vectors.append(vector)
                vector.vector_pools.append(vector_pool)
            
        VECTOR_POOLS.append(vector_pool)
        VECTOR_POOLS_DICT[vector_pool.id] = vector_pool
        
    injections_df = read_injections_df(database_file_paths["Injections"])
    
    for injection_row in injections_df.iterrows():
        
        injection = Injection()
        injection.id = injection_row[0]
        animal = ANIMALS_DICT[injection_row[1]["Animal"]]
        animal.injections.append(injection)
        injection.animal = animal
        quantity = injection_row[1]["Quantity (vgs)"]
        if quantity:
            injection.quantity = float(quantity)
        vector_pool = VECTOR_POOLS_DICT[injection_row[1]["Vector Pool"]]
        vector_pool.injections.append(injection)
        injection.vector_pool = vector_pool
            
        INJECTIONS.append(injection)
        INJECTIONS_DICT[injection.id] = injection
        
    dissociation_runs_df = read_dissociation_runs_df(database_file_paths["Dissociation Runs"])
    
    for dissociation_run_row in dissociation_runs_df.iterrows():
        
        dissociation_run = Dissociation_Run()
        dissociation_run.name = dissociation_run_row[0]
        
        if dissociation_run_row[1]["Samples Processed"] is not None:
            for tissue_sample_name in dissociation_run_row[1]["Samples Processed"].split(","):
                tissue_sample = TISSUE_SAMPLES_DICT[tissue_sample_name]
                dissociation_run.tissue_samples.append(tissue_sample)
                tissue_sample.dissociation_run = dissociation_run
                
        dissociation_run.protocol_version = float(dissociation_run_row[1]["10x Version"])
            
        DISSOCIATION_RUNS.append(dissociation_run)
        DISSOCIATION_RUNS_DICT[dissociation_run.name] = dissociation_run
        
    sequencing_runs_df = read_sequencing_runs_df(database_file_paths["Sequencing Runs"])
    
    for sequencing_run_row in sequencing_runs_df.iterrows():
        
        sequencing_run = Sequencing_Run()
        sequencing_run.name = sequencing_run_row[0]
        
        if sequencing_run_row[1]["Related to Read Sets (Sequencing Run)"] is not None:
            for read_set_name in sequencing_run_row[1]["Related to Read Sets (Sequencing Run)"].split(","):
                read_set = READ_SETS_DICT[read_set_name]
                sequencing_run.read_sets.append(read_set)
                read_set.sequencing_run = sequencing_run
                
        sequencing_run.read_1_length = int(sequencing_run_row[1]["Read 1 Length"])
        sequencing_run.read_2_length = int(sequencing_run_row[1]["Read 2 Length"])
            
        SEQUENCING_RUNS.append(sequencing_run)
        SEQUENCING_RUNS_DICT[sequencing_run.name] = sequencing_run

def read_linked_df(csv_file_path, linked_fields=None):
    
    df = pandas.read_csv(csv_file_path, index_col=0, header=0)
    
    for row in df.iterrows():
        
        if linked_fields is None:
            continue
        
        for linked_field in linked_fields:
            
            if isinstance(row[1][linked_field], str):
                string_values = row[1][linked_field].split(", ")
                df.loc[row[0], linked_field] = ",".join([x for x in string_values])
            else:
                df.loc[row[0], linked_field] = None
            
    return df
    
def read_tissue_samples_df(csv_file_path):
    
    return read_linked_df(csv_file_path, ["Cell Sets", "Dissociation Run", "Animal"])
    
def read_cell_sets_df(csv_file_path):
    
    return read_linked_df(csv_file_path, ["Sequencing Libraries", "Source Tissue"])

def read_read_sets_df(csv_file_path):
    
    return read_linked_df(csv_file_path, ["Libraries", "Sequencing Run"])

def read_sequencing_libraries_df(csv_file_path):
    
    return read_linked_df(csv_file_path, ["Read Sets", "Sequencing Runs", "Cell Set"])

def read_alignments_df(csv_file_path):
    
    return read_linked_df(csv_file_path, ["Reference"])

def read_references_df(csv_file_path):
    
    return read_linked_df(csv_file_path, ["Alignments"])

def read_templates_df(csv_file_path):
    
    return read_linked_df(csv_file_path, ["Vectors"])

def read_vectors_df(csv_file_path):
    
    return read_linked_df(csv_file_path, ["Delivery Vehicle", "Cargo", "Vector Pools"])

def read_animals_df(csv_file_path):
    
    return read_linked_df(csv_file_path, ["Injections", "Samples"])

def read_vector_pools_df(csv_file_path):
    
    return read_linked_df(csv_file_path, ["Vectors", "Injections"])

def read_injections_df(csv_file_path):
    
    return read_linked_df(csv_file_path, ["Animal", "Vector Pool"])

def read_viruses_df(csv_file_path):
    
    return read_linked_df(csv_file_path, ["Vectors"])

def read_dissociation_runs_df(csv_file_path):
    
    return read_linked_df(csv_file_path, ["Samples Processed"])

def read_sequencing_runs_df(csv_file_path):
    
    return read_linked_df(csv_file_path, ["Related to Read Sets (Sequencing Run)"])

def get_data_path():
    
    return DATA_PATH

load_database()

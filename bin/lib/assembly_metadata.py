# assembly metadata

import os
import pandas
import argparse

FLAGS = None
UNPARSED = None


class AssemblyMetadata:

    def __init__(self, *,assembly_folder):
        self.num_assembly = None
        self.sequence_length=None
        self.assembly_length_by_assembly_id = self.load_assembly_length(assembly_folder=assembly_folder)
        self.assembly_path_by_assembly_id = self.load_assembly_path(assembly_folder=assembly_folder)
        self.tax_id_by_assembly_id, self.assembly_id_by_species_tax_id = self.load_assembly_tax_id(assembly_folder=assembly_folder)
        self.assembly_id_by_sequence_id, self.sequence_id_by_assembly_id = self.load_sequence_summary(assembly_folder=assembly_folder)
        self.human_noise_bed_path = self.load_human_noise_bed_path(assembly_folder=assembly_folder)


    def get_num_assembly(self):

        return self.num_assembly

    def get_assembly_length(self, *, assembly_list, how='inner'):

        return assembly_list.merge(
                                   right=self.assembly_length_by_assembly_id,
                                   how=how,
                                   left_on='assembly_id',
                                   right_index = True,
                                   suffixes=['', '_y'],
                                   validate='m:1'
                                  )

    def get_assembly_path(self, *, assembly_list, how='inner'):

        return assembly_list.merge(
                                   right=self.assembly_path_by_assembly_id,
                                   how=how,
                                   left_on='assembly_id',
                                   right_index = True,
                                   suffixes=['', '_y'],
                                   validate='m:1'
                                  )

    def get_tax_id(self, *, assembly_list, how='inner'):

        return assembly_list.merge(
                                   right=self.tax_id_by_assembly_id,
                                   how=how,
                                   left_on='assembly_id',
                                   right_index = True,
                                   suffixes=['', '_y'],
                                   validate='m:1'
                                  )

    def get_sequence_tax_id(self, *, assembly_list, how='inner'):

        assembly_tax_id = self.get_tax_id(assembly_list=assembly_list, how=how)

        return assembly_tax_id.merge(
                                     right=self.sequence_id_by_assembly_id,
                                     how=how,
                                     left_on='assembly_id',
                                     right_index = True,
                                     suffixes=['', '_y'],
                                     validate='m:m'
                                    )

    def get_assembly_id(self, *, species_tax_id_list, how='inner'):

        return species_tax_id_list.merge(
                                         right=self.assembly_id_by_species_tax_id,
                                         how=how,
                                         left_on='species_tax_id',
                                         right_index = True,
                                         suffixes=['', '_y'],
                                         validate='m:m'
                                        )

    def get_sequence_length(self, *, sequence_list, how='inner'):

        return sequence_list.merge(
                                   right=self.assembly_id_by_sequence_id,
                                   how=how,
                                   left_on='sequence_id',
                                   right_index = True,
                                   suffixes=['', '_y'],
                                   validate='m:1'
                                  )

    def get_human_noise_bed_path(self, *, assembly_list, how='inner'):

        return assembly_list.merge(
                                   right=self.human_noise_bed_path,
                                   how=how,
                                   left_on='assembly_id',
                                   right_index = True,
                                   suffixes=['', '_y'],
                                   validate='m:1',
                                  )


    def load_assembly_length(self, *, assembly_folder):

        assembly_length_col = [
                               'assembly_id', 
                               'assembly_length',
                              ];
        assembly_length_col_type = {
                                    'assembly_id': str,
                                    'assembly_length': int,
                                   }

        assembly_length = pandas.read_csv(os.path.join(assembly_folder, 'assembly_length'),
                                          sep='\t',
                                          names=assembly_length_col,
                                          usecols=assembly_length_col,
                                          dtype=assembly_length_col_type,
                                          header=None,
                                          memory_map=True,
                                         )
        
        if self.num_assembly is None:
            self.num_assembly = assembly_length.shape[0]
        else:
            if assembly_length.shape[0] != self.num_assembly:
                print('Number of records in assembly_length is not consistent', file=os.sys.stderr)
        
        return assembly_length.set_index('assembly_id')


    def load_assembly_path(self, *, assembly_folder):

        assembly_path_col = [
                             'assembly_id', 
                             'path',
                            ];
        assembly_path_col_type = {
                                  'assembly_id': str,
                                  'path': str,
                                 }

        assembly_path = pandas.read_csv(os.path.join(assembly_folder, 'assembly_path'),
                                        sep='\t',
                                        names=assembly_path_col,
                                        usecols=assembly_path_col,
                                        dtype=assembly_path_col_type,
                                        header=None,
                                        memory_map=True,
                                       )

        if self.num_assembly is None:
            self.num_assembly = assembly_path.shape[0]
        else:
            if assembly_path.shape[0] != self.num_assembly:
                print('Number of records in assembly_path is not consistent', file=os.sys.stderr)
        
        return assembly_path.set_index('assembly_id')


    def load_assembly_tax_id(self, *, assembly_folder):

        assembly_tax_id_col = [
                               'assembly_id', 
                               'tax_id',
                               'species_tax_id',
                               'genus_tax_id',
                               'genus_height',
                              ];
        assembly_tax_id_col_type = {
                                    'assembly_id': str,
                                    'tax_id': int,
                                    'species_tax_id': int,
                                    'genus_tax_id': int,
                                    'genus_height': int,
                                   }

        assembly_tax_id = pandas.read_csv(os.path.join(assembly_folder, 'assembly_tax_id'),
                                          sep='\t',
                                          names=assembly_tax_id_col,
                                          usecols=assembly_tax_id_col,
                                          dtype=assembly_tax_id_col_type,
                                          header=None,
                                          memory_map=True,
                                         )
        
        if self.num_assembly is None:
            self.num_assembly = assembly_tax_id.shape[0]
        else:
            if assembly_tax_id.shape[0] != self.num_assembly:
                print('Number of records in assembly_tax_id is not consistent', file=os.sys.stderr)

        return assembly_tax_id.set_index('assembly_id'), assembly_tax_id.set_index('species_tax_id')


    def load_sequence_summary(self, *, assembly_folder):

        sequence_summary_col = [
                                'sequence_id',
                                'sequence_length',
                                'assembly_id',
                               ];
        sequence_summary_col_type = {
                                     'sequence_id': str,
                                     'sequence_length': int,
                                     'assembly_id': str,
                                    }

        sequence_summary = pandas.read_csv(os.path.join(assembly_folder, 'sequence_summary'),
                                           sep='\t',
                                           names=sequence_summary_col,
                                           usecols=sequence_summary_col,
                                           dtype=sequence_summary_col_type,
                                           header=None,
                                           memory_map=True,
                                          )

        return sequence_summary.set_index('sequence_id'), sequence_summary.set_index('assembly_id')


    def load_human_noise_bed_path(self, *, assembly_folder):

        human_noise_bed_path_col = [
                                    'assembly_id', 
                                    'path',
                                    #'span_bp',
                                   ];
        human_noise_bed_path_col_type = {
                                         'assembly_id': str,
                                         'path': str,
                                         #'span_bp': int,
                                        }

        human_noise_bed_path = pandas.read_csv(os.path.join(assembly_folder, 'human_noise_bed_path'),
                                               sep='\t',
                                               names=human_noise_bed_path_col,
                                               usecols=human_noise_bed_path_col,
                                               dtype=human_noise_bed_path_col_type,
                                               header=None,
                                               memory_map=True,
                                              )
        
        return human_noise_bed_path.set_index('assembly_id')


    def data_integrity_check(self):

        data_integrity_check_passed = True

        if self.assembly_length_by_assembly_id.shape[0] != self.assembly_path_by_assembly_id.shape[0] or self.assembly_length_by_assembly_id.shape[0] != self.tax_id_by_assembly_id.shape[0]:
            data_integrity_check_passed = False
            print('AssemblyMetadata - Number of assembly_id not match', file=os.sys.stderr, flush=True)


        assembly = self.assembly_length_by_assembly_id.merge(
                                                             right=self.assembly_path_by_assembly_id,
                                                             how='inner', 
                                                             left_index=True, 
                                                             right_index = True,
                                                             suffixes=['', '_y'],
                                                             validate='1:1'
                                                            )
        if assembly.shape[0] != self.assembly_length_by_assembly_id.shape[0]:
            data_integrity_check_passed = False
            print('Some assembly_id cannot be matched - assembly_length vs assembly_path', file=os.sys.stderr, flush=True)

        assembly = self.assembly_length_by_assembly_id.merge(
                                                             right=self.tax_id_by_assembly_id,
                                                             how='inner', 
                                                             left_index=True, 
                                                             right_index = True,
                                                             suffixes=['', '_y'],
                                                             validate='1:1'
                                                            )
        if assembly.shape[0] != self.assembly_length_by_assembly_id.shape[0]:
            data_integrity_check_passed = False
            print('Some assembly_id cannot be matched - assembly_length vs assembly_tax_id', file=os.sys.stderr, flush=True)

        assembly = self.assembly_path_by_assembly_id.merge(
                                                           right=self.tax_id_by_assembly_id,
                                                           how='inner', 
                                                           left_index=True, 
                                                           right_index = True,
                                                           suffixes=['', '_y'],
                                                           validate='1:1'
                                                          )
        if assembly.shape[0] != self.assembly_path_by_assembly_id.shape[0]:
            data_integrity_check_passed = False
            print('Some assembly_id cannot be matched - assembly_path vs assembly_tax_id', file=os.sys.stderr, flush=True)


        assembly_list = self.sequence_length[['assembly_id']].drop_duplicates().set_index(['assembly_id'])

        assembly = self.assembly_length_by_assembly_id.merge(
                                                             right=assembly_list,
                                                             how='inner', 
                                                             left_index=True, 
                                                             right_index = True,
                                                             suffixes=['', '_y'],
                                                             validate='1:1'
                                                            )
        if assembly.shape[0] != assembly_list.shape[0]:
            data_integrity_check_passed = False
            print('Some assembly_id cannot be matched - sequence_summary vs assembly_length', file=os.sys.stderr, flush=True)

        
        if data_integrity_check_passed == True:
            print('Data integrity check passed.')
            print('Number of assembly = {num_assembly}'.format(num_assembly=self.assembly_length_by_assembly_id.shape[0]))
            print('Number of sequence = {num_sequence}'.format(num_sequence=self.sequence_length.shape[0]))
        else:
            print('Data integrity check failed.')
            print('assembly_length : {assembly_length}'.format(assembly_length=self.assembly_length_by_assembly_id.shape[0]))
            print('assembly_path   : {assembly_path}'.format(assembly_path=self.assembly_path_by_assembly_id.shape[0]))
            print('assembly_tax_id : {assembly_tax_id}'.format(assembly_tax_id=self.tax_id_by_assembly_id.shape[0]))
            print('sequence_length : {sequence_length}'.format(sequence_length=self.sequence_length.shape[0]))


def main():

    assembly_metadata = AssemblyMetadata(assembly_folder=FLAGS.assemblyFolder)
    assembly_metadata.data_integrity_check()



if __name__ == '__main__':
    nano_dir=os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    parser = argparse.ArgumentParser(description='assembly_metadata')

    parser.add_argument('--assemblyFolder', help='Assembly folder', default=nano_dir+'/genomes')

    FLAGS, UNPARSED = parser.parse_known_args()

    main()

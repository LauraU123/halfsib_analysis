import pandas as pd
import argparse

#f = 'example2/filtered_SNP_map.txt'

#outfile = 'example2/plink.map'

def SNP_map(f, outfile):

	print('Reading in file')
	df = pd.read_csv(f, sep='\t')
	print('Read in file')
	
	df['Name'] = list( map(lambda x: x.upper(), df['Name']) )
	print(df)
	print('Making SNP Name uppercase')
	
	df.insert( 2,'Distance', '0')
	#print(df)
	print('Insert column Distance filled with O')
	
	dfsub = df[ [ 'Chromosome', 'Name', 'Distance', 'Position' ] ]
	print('re-ordered columns')

	# if X chromosome is there
	#dfsub.drop(df[df['Chromosome'] == 'X'].index, inplace = True)
	#dfsub.drop(df[df['Chromosome'] == 'MT'].index, inplace = True)
	#dfsub.drop(df[df['Chromosome'] == 'Y'].index, inplace = True)
	#dfsub.drop(df[df['Chromosome'] == 0].index, inplace = True)
	
	dfsub.to_csv(outfile, index=False, header=False, sep='\t') # need header as false to NOT output number indices
	print('wrote to file: %s' % (outfile))


#f = 'example2/Sample_Map.txt'

#outfile = 'example2/plink.fam'

def Sample_map(f, outfile):

	print('Reading in file')
	df = pd.read_csv(f, sep='\t')
	print('Read in file')
	
	df.insert(3, 'Zero', '0')
	print('Insert column Distance filled with O')
	
	dfsub = df[ [ 'Name', 'ID', 'Zero', 'Zero', 'Zero', 'Zero'] ]
	print('re-ordered columns')
	dfsub['Name'] = dfsub['Name'].astype(str).str.replace(' ','_')
	dfsub['ID'] = dfsub['ID'].astype(str).str.replace(' ','_')
	dfsub.to_csv(outfile, index=False, header=False, sep='\t') # need header as false to NOT output number indices
	print('wrote to file: %s' % (outfile))


#f = 'example2/Univ_of_Bern_Drogemuller_BOV770V01_20191024_FinalReport.txt'

#outfile = 'example2/plink.lgen'

def lgen(f, outfile):
	
	print('Reading in file')
	df = pd.read_csv(f, skiprows=[0,1,2,3,4,5,6,7,8], sep='\t')
	print('Read in file and skipped first 10 rows (i.e. the header)')
	#df.drop(df.index[[0,1,2,3,4,5,6,7,8,9]])
	
	df['Allele1 - Forward'] = df['Allele1 - Forward'].str.replace('-','0')
	df['Allele2 - Forward'] = df['Allele2 - Forward'].str.replace('-','0')
	df['Sample ID'] = df['Sample ID'].astype(str).str.replace(' ','_')
	
	print('replaced - with 0 in genotypes Allele1 - Forward and Allele2 - Forward')
	
	df['SNP Name'] = list( map(lambda x: x.upper(), df['SNP Name']) )
	print('Making SNP Name uppercase')
	
	dfsub = df[ [ 'Sample ID', 'Sample ID', 'SNP Name', 'Allele1 - Forward', 'Allele2 - Forward' ] ]
	print('re-ordered columns')
	
	dfsub.to_csv(outfile, index=False, header=False, sep=' ') # need header as false to NOT output number indices
	print('wrote to file: %s' % (outfile))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Generate Plink Input files",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--SNP_in", required=True, help=".map file")
    parser.add_argument("--SNP_out", required=True, help=".map file")
    parser.add_argument("--sample_in", required=True, help=".map file")
    parser.add_argument("--sample_out", required=True, help=".csv file")
    parser.add_argument("--experiment_in", required=True, help=".map file")
    parser.add_argument("--experiment_out", required=True, help=".map file")
    args = parser.parse_args()
    
	
    SNP_map(args.SNP_in, args.SNP_out)
    Sample_map(args.sample_in, args.sample_out)
    lgen(args.experience_in, args.experience_out)
#!/usr/bin/env python
# coding: utf-8

#------------------------------
# Import the needed libraries 
#------------------------------
import pandas as pd
import numpy as np
import csv, os, sys
import logging
import argparse
from sklearn.model_selection import train_test_split
# Import DecisionTreeClassifier
from sklearn.tree import DecisionTreeRegressor
# Import AdaBoostClassifier
from sklearn.ensemble import AdaBoostRegressor
from sklearn.model_selection import GridSearchCV
import xgboost as xgb


# max depth 3 
# nestimators 20,000
# cv ten times


logging.basicConfig(level=logging.INFO, 
                    format='%(asctime)s : %(levelname)s : %(message)s',
                    datefmt='%H:%M:%S')
logger = logging.getLogger(__name__)


class XGBoost:
    '''
    Define the Adaboost object 
    '''
    def __init__(self, args):
        '''
        Create properties for the new object and assign values to them
            Parse the arguments given to the script 
            should include: The VCF file of interest 
                            HapMap file destination 
                            Output file location 
        '''
        
        #***********************************
        # create and assign the attributes 
        #***********************************
        self.vcf_path = args.vcf 
        self.truth_path = args.truth_vcf
        self.output_path = args.output_dir
        self.features = args.features.split(",")

        if not os.path.exists(self.output_path):
            os.makedirs(self.output_path)

        #------------------------------------------
        # Read in the VCF file as a pandas dataframe and store it in the object
        #------------------------------------------
        self.vcf = pd.read_csv(self.vcf_path,
                         sep='\t', low_memory=False, comment='#', header =None,
                         names=["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "ENCODING"])

        #------------------------------------------
        # Read in the HapMap information and store it in the object 
        #------------------------------------------
        self.truth = pd.read_csv(self.truth_path, 
                        sep='\t', header=None, comment= "#", low_memory=False,
                       names=["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "ENCODING"])



    #********************************************
    # Additional methods for the Adaboost Object 
    #********************************************
    
    def configure_vcf(self):
        '''
        Function to read in and configure the VCF files
            Configure the VCF and separate the lines into a dictionary by chromosome 
            Return the VCF header and dictionary
        '''
        vcf_header = []
        vcf_dic = {}
        vcf = open( self.vcf_path, "r" )
        for line in csv.reader(vcf, delimiter = "\t"):
            # check to see if header, is so append the line to header 
            if line[0][0] == "#":
                vcf_header.append(line)
                continue
            # else append the line to the dictionary based on the chromosome
            chromosome = line[0]
            position = line[1]
            new_key = str(chromosome) + ":" + str(position)
            vcf_dic[new_key]= line
        vcf.close()

        # Add the vcf header and body to the new object 
        self.vcf_header = vcf_header
        self.vcf_body = vcf_dic
        return( self )

    def abaBoost(self, X, y, X_train, X_test, y_train, y_test):
        '''
        Function to create the Gradient boosting model.
        Takes as input the feature data, data labels, and the data split into training and test data sets. Outputs the labels given to 
        the data predicted by the model.
            Inputs
                X: Feature data
                y: labels given to the data
                X_train: the feature data assigned to the training data
                X_test: the feature. data assigned to the test data set 
                y_train: labels that are assigned to the traing data set 
                y_test: labels that are assigned to the test data set 
            output
                y_predicted: The labels predicted by the model 
        '''
        logger.info("Running Gradient boosting regression.")
        '''

        tree_depth = interaction_depth + 1

        max_depth:         Maximum depth of the tree
        min_samples_split: the minimum number of samples required to split an internal node (default = 2)
        min_samples_leaf:  Minimum number of samples required to be a leafe node. May effect the smoothing of the model, especially with regression (def=1)
        min_weight_fraction_leaf:
        '''

        #-------------------------------------
        # Run Gradient boosting regression 
        #-------------------------------------
        xg_regression_model = xgb.XGBRegressor(
                              gamma = 1,
                              max_depth = 1,
                              min_child_weight = 27,
                              n_estimators=20000
                              )


        # Fit training data to the model
        xg_regression_model.fit(X_train, y_train)
        # Get the predicted values 
        y_predicted = xg_regression_model.predict(X)

        self.data['xgb_score'] = y_predicted

        return self 


    def preprocessingData(self):
        '''
        Load in the data VCF then convert info column into pandas dataframe 
        '''

        logger.info("Processing the input data.")
        
        ## subset the data to get the get 'Chr', 'Pos','REF','ALT'
        df_vcf_subset = self.vcf[['CHROM', 'POS','REF','ALT']]


        #---------------------------------------------------------------------
        # Load in the data VCF and convert info column into pandas dataframe 
        #---------------------------------------------------------------------
        ## Read in the data header of the vcf file to get info column ID's
        ## Seperate the info column of the vcf for each variant 
        ## create a pandas dataframe for the information column 


        # lists to hold the ID's and the values 
        info_id = []
        all_info = []
        # Get the ID's from the header
        for i in self.vcf_header:
            if i[0][0:11] == "##INFO=<ID=":
                info_id.append(str(i[0].split(",")[0][11:])) 
        # print(info_id)

        # Iterate through each variant
        for i in self.vcf_body:
            info_num = [None]*len(info_id)
            ## split the info section
            info = self.vcf_body[i][7].split(";")
            if "" in info:
                info.remove("")
            for value in info:
                ## pull out the ID and value 'ID=Value'
                temp = value.split("=")
                ## If the ID has no value (given by IndexError), make binary indicator, 1=pressent 
                try: 
                    info_num[info_id.index(temp[0])] = temp[1]
                except IndexError:
                    info_num[info_id.index(temp[0])] = 1
            all_info.append(info_num)
        df_info = pd.DataFrame(data = all_info)
        df_info.columns = info_id
        print(df_info.head())

        #------------------------------------------------------------------
        # Combine the new INFO DataFrame (df_info) with the VCF (df_vcf)
        #------------------------------------------------------------------
        ## Combine the VCF's first 5 columns (CHROM, POS, ID, ALT, REF), to the newly created INFO dataframe 
        DF = pd.concat([self.vcf.iloc[:,0:5], df_info], axis=1)
        ## Create the unique index column "CHROM:POS"
        DF['IND'] = DF['CHROM'] + ':' + DF['POS'].astype(str)
        ## Creeate the SNP "REF:ALT"
        DF['SNP'] = DF['REF'] + ':' + DF['ALT']
        ## Set the GT from the VCF, (ex. "1/1")
        DF['GT'] = self.vcf.iloc[:,-1].apply(lambda x: x.split(':')[0])
        ## Set the index for the DataFrame
        DF = DF.set_index('IND')

        #-----------------------------------------
        # Remove indels - only SNPs are retained
        #-----------------------------------------
        ## index to variants with ALT and REF bases equal to one
        DF = DF[DF['ALT'].str.len().eq(1)]
        DF = DF[DF['REF'].str.len().eq(1)]

        print(DF.head())

        #---------------------------------------------
        #Subset the DataFrame to select info features 
        #---------------------------------------------
        DF_subset = DF[['AC', 'AF', 'AN', 'BaseQRankSum', 'DP', 'Entropy', 'FS', 
               'Homopolymer', 'MLEAC', 'MLEAF', 'MMF', 'MQ', 'MQRankSum',
               'QD', 'RNAEDIT', 'RPT', 'ReadPosRankSum',  'SOR', 'SPLICEADJ',
               'TDM', 'TMMR', 'TPR', 'VAF','VPR', 'SNP','GT', 'DJ', 'ED', 'RS']]

        print(DF_subset.head())

        #-------------------------------------------
        # Edit the info outputs 
        #-------------------------------------------
        ## Set all None's to 0's
        DF_subset = DF_subset.fillna(0)
        ## Set values from None to 0, and if the value is not 0 then set to 1
        DF_subset['RNAEDIT'][DF_subset['RNAEDIT'] != 0] = 1
        DF_subset['RPT'][DF_subset['RPT'] != 0] = 1

        ## Set the GT's to 0 if "0/1" and 1 if "1/1"
        DF_subset = pd.get_dummies(DF_subset, columns=['GT'], prefix='GT')

        ## pd.get_dummies creates a binary matrix where values in a list/vector are similar or not 
        DF_subset = pd.get_dummies(DF_subset, columns=['SNP'], prefix='SNP')

        # Set the splice junction distance to 1 if the value is greater then 0
        # DF_subset[DF_subset['SPLICEADJ'].astype(int)>0] = 1

        print(DF_subset.head())

        #-------------------------------------------
        # Subset and set the values to floats
        #-------------------------------------------
        # DF_Features = DF_subset[['AC', 'AF', 'AN', 'BaseQRankSum', 'DP', 'Entropy', 'FS', 'Homopolymer',
        #    'MLEAC', 'MLEAF', 'MMF', 'MQ', 'MQRankSum', 'QD', 'RNAEDIT', 'RPT',
        #    'ReadPosRankSum', 'SOR', 'SPLICEADJ', 'TDM', 'TMMR', 'TPR', 'VAF',
        #    'VPR', 'DJ', 'ED', 'GT_0/1', 'GT_1/1', 'SNP_A:C', 'SNP_A:G', 'SNP_A:T',
        #    'SNP_C:A', 'SNP_C:G', 'SNP_C:T', 'SNP_G:A', 'SNP_G:C', 'SNP_G:T',
        #    'SNP_T:A', 'SNP_T:C', 'SNP_T:G', 'RS']]
        # DF_Features = DF_subset[["QD","ReadPosRankSum", "FS", "VPR", "VAF", "SPLICEADJ", "RPT", "Homopolymer", "Entropy", "RNAEDIT"]]
        print(self.features)
        DF_Features = DF_subset[self.features]

        # DF_Features = DF_Features[DF_Features['RPT']==0.0]
        # DF_Features = DF_Features.drop(columns='RPT')
        DF_Features = DF_Features.astype(float)

        DF_Features.to_csv('DF_Features.vcf', sep='\t')

        self.DF_Features = DF_Features

        return self


    def processTruthSet(self):
        #-------------------------------------------
        ## Load in the HAPMAP data set 
        #-------------------------------------------
        logger.info("Processing the truth set file.")
        df_truth = self.truth

        # # if only the ids are given
        # if all(df_truth[["POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "ENCODING"]].isnull()):
        #     # make sure the column are strings 
        #     df_truth['IND'] = df_truth['CHROM'].map(str)
            
        # else:
        # Create the unique index column for each variant CHR:POS
        df_truth['IND'] = df_truth['CHROM']+':'+ df_truth['POS'].astype(str)
        #--------------------------------------------
        # Identify if the variants are in Truth Set
        #--------------------------------------------
        # Create a column of 1's 
        df_truth['In_TRUTHSET'] = np.ones(df_truth.shape[0])
        # Reset the index for the dataframe 
        df_truth = df_truth.set_index('IND')
        df_truth.index.name = 'IND'
        print(df_truth)

        # Merge the subsetted VCF with HAPMAP 
        df_train = self.DF_Features.merge( df_truth['In_TRUTHSET'], 
                                        how='left', left_index=True, right_index=True)
        df_train['In_TRUTHSET'] = df_train['In_TRUTHSET'].fillna(0)
        print(df_train.head())
        df_train.to_csv('DF_Train.vcf', sep='\t')

        self.data = df_train
        
        return self

def ecdf(data):
    '''
    Function to calculate the Empirical Cumulative Distribution Function (ECDF)
    '''
    """ Compute ECDF """
    # x = np.sort(data)
    n = data.shape[0]
    y = np.arange(1, n+1) / n
    print(data)
    data["ecdf"] = y
    
    return data 

def splitData(obj):
    #---------------------------------------
    # split the dataset
    #---------------------------------------
    logger.info("Splitting the data into test data.")
    # X is the feature data 
    X = obj.data[obj.data.columns[:-1]]
    # y is the label data
    y = obj.data['In_TRUTHSET']
    ## Choose test size or do Cross validation
    X_train, X_test, y_train, y_test = train_test_split(X, y, 
                                                        test_size=0.5, 
                                                        stratify= y, 
                                                        random_state=1)
    return X, y, X_train, X_test, y_train, y_test




def predictedOutput(obj):
    '''
    Get ECCDF values for the xgbboost scores
    
    Filter the variants to only include common variants set by the DBSNP value RS
    Get the ECDF values for each of these variants using their given xgbboost scores 
    Remove the variants with an ECDF value bellow 0.05
        Inputs
            obj    : The ADABoost class object 
            y_pred : The predicted xgbboost values 
            X      : Feature dataset 
    '''

    # get variants that are common variants 
    RS_subset = obj.data[ obj.data['In_TRUTHSET'] == 1 ]
    # RS_subset = obj.data

    # sort the variants by their xgbboost score 
    RS_subset = RS_subset.sort_values(by = "xgb_score")

    logger.info("Generating ECDF values")

    # Generate the ECDF values using the xgbboost scores 
    ecdf_df = ecdf(RS_subset)
    print(ecdf_df.head())
    ecdf_df.to_csv("ecdf_data.txt", sep='\t', index=True)


    
    #-----------------
    # Plot the ECDF
    #-----------------
    logger.info("plotting the ECDF")
    import matplotlib.pyplot as plt
    
    plt.scatter(x = ecdf_df['xgb_score'], 
                y = ecdf_df['ecdf']);
    plt.xlabel('x', fontsize=16)
    plt.ylabel('y', fontsize=16)
    plt.savefig('ecdf_plot.png')
    
    # full.to_csv("ecdf_data.txt", sep='\t', index=False)

    # real_snps = set( ecdf_df[ecdf_df['ecdf'] > 0.05].index )
    real_neg_snps = ecdf_df[ ecdf_df['ecdf'] <= 0.05 ]

    # find the max value for the xgbscore among all the variants that the ECDF value falls 0.05 and bellow 
    max_score = max(real_neg_snps['xgb_score'])
    print("MAX VAL", max_score)

    # pull out the index for the variants that have sgb scores that fall above the 0.05 value 
    real_snps = set( obj.data[ obj.data['xgb_score'] > max_score ].index )
    print("length:",len(real_snps))

    # Get the ID'S CHROM:POS
    obj.vcf['chr:pos'] = obj.vcf['CHROM']+':'+obj.vcf['POS'].astype(str)
    
    # Pull out the identified variants 
    df_bm = obj.vcf[ obj.vcf['chr:pos'].isin(real_snps) ]
    
    # Remove the last column "in_hapmap"
    df_bm = df_bm.iloc[:, :-1]

    ## Save vcf without header
    out_file = os.path.join(obj.output_path, 'xgbboost.vcf')

    # add a header to the vcf 
    header = []
    for i in obj.vcf_header:
        header.append("\t".join(i))
    header = "\n".join(header) + "\n"
    add_header = open(out_file,"w")
    add_header.write(header)
    add_header.close()

    df_bm.to_csv(out_file, sep='\t', index=False, header=False, mode = "a")


def main():
    ###########################################
    # Gather the arguments passed to the SNPiR script in command line 
    ###########################################

    ## Input Arguments
    # Description 
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
        description = "Run XGBoost.\n")
    # Mandatory arguments 
    parser.add_argument('--vcf', required = True,
        help="input vcf.")
    parser.add_argument('--truth_vcf', required = True,
        help="input truth set data.")
    parser.add_argument("--output_dir", required=True, 
        help="output directory")
    parser.add_argument("--features", 
                        required = False, 
                        type     = str,
                        help     = "Which features to use when boosting.",
                        default = "QD,ReadPosRankSum,FS,VPR,VAF,VMMF,SPLICEADJ,RPT,Homopolymer,Entropy,RNAEDIT")
    
    # Parse the given arguments 
    args = parser.parse_args()
    
    logger.info("Initializing XGBoost class object!")
    xgb_object = XGBoost(args)

    xgb_object = xgb_object.configure_vcf()
    
    xgb_object = xgb_object.preprocessingData()

    xgb_object = xgb_object.processTruthSet()

    

    #---------------------------------------
    # split the dataset
    #---------------------------------------
    X, y, X_train, X_test, y_train, y_test = splitData(xgb_object)

    print(X_train.shape, X_test.shape, y_train.shape, y_test.shape)
    print(y_train.value_counts(),'\n', y_test.value_counts())

    # xgb_object = abaBoost(xgb_object, X, y, X_train, X_test, y_train, y_test)
    xgb_object = xgb_object.abaBoost(X, y, X_train, X_test, y_train, y_test)


    predictedOutput( obj = xgb_object )


   # sys.exit(0)

if __name__ == "__main__":

    main()

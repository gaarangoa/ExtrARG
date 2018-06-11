import numpy as np
import pandas as pd
from sklearn.feature_selection import SelectFromModel
from sklearn.ensemble import ExtraTreesClassifier as ETC
from sklearn.cross_validation import cross_val_score
from bayes_opt import BayesianOptimization
import click

class ExtraARG():
    def __init__(self, input_file, output_file, disc, min_reads):
        self.input_file = input_file
        self.output_file = output_file
        self.disc = disc
        self.min_reads = min_reads
        self.X = []
        self.y = []
        self.z = []

    def load_data(self):
        Raw_data2=pd.read_excel(self.input_file,sheet_name='Raw')
        Raw_data=pd.read_excel(self.input_file,sheet_name='Scale')
        Map=pd.read_excel(self.input_file,sheet_name='Map')
        Raw_data.index=Raw_data.pop('samples')

        listR=[]
        listD=[]
        for i in range(0, len(Raw_data.index)):
            counter=0
            for j in range(0,len(Raw_data.columns)):
                if (Raw_data.iloc[i][j]<self.min_reads):
                    counter=counter+1;
            if (counter!=len(Raw_data.columns)):
                listR.append(Raw_data.index[i])
            else:
                listD.append(Raw_data.index[i])

        Raw_data2.index=Raw_data2.pop('samples')
        Raw_data2=Raw_data2.loc[listR]
        Raw_data2=Raw_data2.transpose()
        Map.index=Map.pop('Name')
        Map=Map.loc[Raw_data2.index]
        self.z=Map.Samples
        self.y=Map.Label
        self.X = Raw_data2

    def etc_ccv(self, n_estimators, min_samples_split, max_features):
        val = cross_val_score(
            ETC(n_estimators=int(n_estimators),
                min_samples_split=int(min_samples_split),
                max_features=min(max_features, 0.999),
                random_state=2
            ),
            self.X, self.y, 'accuracy', cv=2
        ).mean()
        return val

    def optimize(self):
        self.gp_params = {"alpha": 1e-5}
        self.etc_0 = BayesianOptimization(
            self.etc_ccv,
            {
                'n_estimators': (10, 1000),
                'min_samples_split': (2, 25),
                'max_features': (0.1, 0.999)
            }
        )
        self.etc_0.maximize(n_iter=10, **self.gp_params)

        print('performing extra trees classifier ...')
        self.forest = ETC(
            n_estimators=int(self.etc_0.res['max']['max_params']['n_estimators']),
            min_samples_split=int(self.etc_0.res['max']['max_params']['min_samples_split']),
            max_features=min(self.etc_0.res['max']['max_params']['max_features'], 0.999),
            random_state=0
        )
        self.forest.fit(self.X, self.y)

    def discriminate(self):

        self.importances = pd.DataFrame({'feature':self.X.columns,'importance':np.round(self.forest.feature_importances_,10)})
        self.importances = self.importances.sort_values('importance',ascending=False).set_index('feature')
        # print(importances)
        # importances.plot.bar()

        print('selecting discriminatory ARGs ...')
        #selecting number of ARGs based on the model
        self.model=SelectFromModel(self.forest,prefit=True)
        self.X_new=self.model.transform(self.X)
        # print(X_new.shape)
        self.list_Country_Wise=self.importances[0:len(self.X_new.T)]

        #selecting number of ARGs based on user defined number
        self.z=self.importances.index[0:self.disc]
        self.Y = self.X[self.z]
        print('saving into '+self.output_file+'.xlsx file ...')
        #writing in new excel file
        writer = pd.ExcelWriter(self.output_file+'.xlsx')
        self.X.to_excel(writer,'cleaned_raw_input')
        self.Y.to_excel(writer,'discriminatory_args')
        writer.save()

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

@click.command(context_settings = CONTEXT_SETTINGS)
@click.option('--input-file', default='', help='input excel file')
@click.option('--output-file', help='output file where to store the results')
@click.option('--disc', default=50, help='top N discriminative ARGs (default 50)')
@click.option('--min-reads', default=1, help='minimum number of reads on each ARG (default 1)')
@click.option('--optimize', default=False, help='minimum number of reads on each ARG (default 1)')
def process(input_file='', output_file='', disc='', min_reads='', optimize=False):
    """
    This program subtract the top N (50 default) discriminatory antibiotic resistance genes from a set of metagenomics samples.
    Hyperparameters of the supervised machine learning algorithm (extra tree classifier) are automatically tuned using the bayesian optimization.
    """

    if not input_file:
        print('\nUsage: extrarg --help\n')
        exit()

    extra_arg = ExtraARG(input_file, output_file, disc, min_reads)

    print('loading input datasets ...')
    extra_arg.load_data()

    print('performing optimization ...')
    extra_arg.optimize()

    print('building model with optimized parameters ...')
    extra_arg.discriminate()


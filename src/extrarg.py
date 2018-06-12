import numpy as np
import pandas as pd
from sklearn.feature_selection import SelectFromModel
from sklearn.ensemble import ExtraTreesClassifier as ETC
from sklearn.cross_validation import cross_val_score
from bayes_opt import BayesianOptimization
import click


class ExtraARG():
    def __init__(self, input_file, output_file, disc, min_reads, epochs, max_importance, min_importance):
        self.input_file = input_file
        self.output_file = output_file
        self.disc = disc
        self.epochs = epochs
        self.min_reads = min_reads
        self.max_importance = max_importance
        self.min_importance = min_importance
        self.X = []
        self.y = []
        self.z = []

    def load_data(self):
        Raw_data2 = pd.read_excel(self.input_file, sheet_name='Raw')
        Raw_data = pd.read_excel(self.input_file, sheet_name='Scale')
        Map = pd.read_excel(self.input_file, sheet_name='Map')
        Raw_data.index = Raw_data.pop('samples')

        listR = []
        listD = []
        for i in range(0, len(Raw_data.index)):
            counter = 0
            for j in range(0, len(Raw_data.columns)):
                if (Raw_data.iloc[i][j] < self.min_reads):
                    counter = counter+1
            if (counter != len(Raw_data.columns)):
                listR.append(Raw_data.index[i])
            else:
                listD.append(Raw_data.index[i])

        Raw_data2.index = Raw_data2.pop('samples')
        Raw_data2 = Raw_data2.loc[listR]
        Raw_data2 = Raw_data2.transpose()
        Map.index = Map.pop('Name')
        Map = Map.loc[Raw_data2.index]
        self.z = Map.Samples
        self.y = Map.Label
        self.X = Raw_data2

    def etc_ccv(self, n_estimators, max_features, top_args):
        _transition_model = ETC(n_estimators=int(n_estimators),
                                max_features=min(max_features, 0.999),
                                random_state=2
                                )
        _transition_model.fit(self.X, self.y)
        _transition_model_select = SelectFromModel(
            _transition_model, prefit=True, threshold=top_args
        )
        _transition_model_x = _transition_model_select.transform(self.X)

        # In this case there are not features to select, therefore, return 0 as accuracy.
        if _transition_model_x.shape[1] == 0:
            return 0

        val = cross_val_score(
            ETC(n_estimators=int(n_estimators),
                max_features=min(max_features, 0.999),
                random_state=2
                ),
            _transition_model_x, self.y, 'accuracy', cv=2
        ).mean()

        return val

    def optimize(self):
        self.gp_params = {"alpha": 1e-5}
        self.etc_0 = BayesianOptimization(
            self.etc_ccv,
            {
                'n_estimators': (100, 1000),
                'max_features': (0.1, 0.5),
                'top_args': (self.min_importance, self.max_importance)
            }
        )
        self.etc_0.maximize(n_iter=self.epochs, **self.gp_params)

        print('performing extra trees classifier ...')
        self.forest = ETC(
            n_estimators=int(
                self.etc_0.res['max']['max_params']['n_estimators']),
            max_features=min(
                self.etc_0.res['max']['max_params']['max_features'], 0.999),
            random_state=0
        )
        self.forest.fit(self.X, self.y)
        self._selected_features_model = SelectFromModel(
            self.forest, prefit=True, threshold=self.etc_0.res['max']['max_params']['top_args']
        )

        self.x_t_selected = self._selected_features_model.transform(self.X)
        # self.forest.fit(self.x_t_selected, self.y)

        self.x_selected = pd.DataFrame(
            data=self.x_t_selected,
            index=self.X.index,
            columns=self.X.columns[self._selected_features_model.get_support()]
        )

        self.importances = pd.DataFrame({
            'Gene': self.x_selected.columns,
            'importance': self.forest.feature_importances_[self._selected_features_model.get_support()]
        })

    def discriminate(self):

        # self.importances = pd.DataFrame(
        #     {
        #         'feature': self.x_selected.columns,
        #         'importance': np.round(self.forest.feature_importances_, 10)
        #     }
        # )
        # self.importances = self.importances.sort_values(
        #     'importance', ascending=False).set_index('feature')

        print('selecting discriminatory ARGs ...')
        # selecting number of ARGs based on the model
        # self.model = SelectFromModel(self.forest, prefit=True)
        # self._reduced_x = self.model.transform(self.X)
        # print("Samples - Features", self.X_new.shape)

        # selecting number of ARGs based on user defined number
        # self.z = self.importances.index
        # self._reduced_x = self.x_selected[self.z]
        # print("Selected ARGs: ", self._reduced_x.shape[1])
        print('saving into '+self.output_file+'.xlsx file ...')
        # writing in new excel file
        writer = pd.ExcelWriter(self.output_file+'.xlsx')
        self.X.to_excel(writer, 'cleaned_raw_input')
        self.x_selected.to_excel(writer, 'discriminatory_args')
        self.importances.to_excel(writer, 'ARGs_importance')
        writer.save()


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('--input-file', default='', help='input excel file')
@click.option('--output-file', help='output file where to store the results')
@click.option('--disc', default=50, help='top N discriminative ARGs (default 50)')
@click.option('--min-reads', default=1, help='minimum number of reads on each ARG (default 1)')
@click.option('--epochs', default=10, help='number of iterations the optimization algorithm run (default 10)')
@click.option('--max-importance', default=0.01, help='maximum importance for search space (default 0.01)')
@click.option('--min-importance', default=1e-5, help='minimum importance for search space (default 1e-6)')
# @click.option('--optimize', default=False, help='minimum number of reads on each ARG (default 1)')
def process(input_file='', output_file='', disc='', min_reads='', epochs=10, max_importance=0.01, min_importance=1e-5):
    """
    This program subtract the top N (50 default) discriminatory antibiotic resistance genes from a set of metagenomics samples.
    Hyperparameters of the supervised machine learning algorithm (extra tree classifier) are automatically tuned using the bayesian optimization.
    """

    if not input_file:
        print('\nUsage: extrarg --help\n')
        exit()

    extra_arg = ExtraARG(input_file, output_file, disc,
                         min_reads, epochs, max_importance, min_importance)

    print('loading input datasets ...')
    extra_arg.load_data()

    print('performing optimization ...')
    extra_arg.optimize()

    print('building model with optimized parameters ...')
    extra_arg.discriminate()

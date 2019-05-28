import numpy as np
import pandas as pd
from sklearn.feature_selection import SelectFromModel
from sklearn.ensemble import ExtraTreesClassifier as ETC
from sklearn.model_selection import cross_validate
from sklearn.metrics import fbeta_score, make_scorer
from sklearn.metrics import confusion_matrix
from bayes_opt import BayesianOptimization
import click
from sklearn.metrics import average_precision_score
from sklearn.metrics import precision_recall_curve
import json

from sklearn.cluster import AffinityPropagation
from sklearn import metrics


class ExtraARG():
    def __init__(self, input_file, output_file, min_reads, epochs, max_importance, min_importance):
        self.input_file = input_file
        self.output_file = output_file
        self.epochs = epochs
        self.min_reads = min_reads
        self.max_importance = max_importance
        self.min_importance = min_importance
        self.X = []
        self.y = []
        self.z = []
        self.cumulative_objective_function = []

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
                if (Raw_data.iloc[i][j] <= self.min_reads):
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
        self._y = Map.Label
        self.map_y = {i: ix for ix, i in enumerate(list(set(self._y)))}
        self.y = np.array([self.map_y[i] for i in self._y])
        self.X = Raw_data2

    def loss_function(self, y_test, y_score):
        cm = confusion_matrix(y_test, y_score)
        precisions = []
        for ix, row in enumerate(cm):
            try:
                precisions.append(float(row[ix]) / sum(row))
            except:
                precisions.append(0)
        precisions = np.array(precisions)

        recalls = []
        for ix, row in enumerate(np.transpose(cm)):
            try:
                recalls.append(float(row[ix]) / sum(row))
            except:
                recalls.append(0)

        recalls = np.array(recalls)
        return 2*precisions.mean()*recalls.mean()/(sum(precisions)+sum(recalls))

    def loss_function_2(self, y_test, y_score):
        # print(y_test, y_score)
        return metrics.adjusted_rand_score(y_test, y_score)

    def etc_ccv(self, n_estimators, top_args):
        _transition_model = ETC(n_estimators=int(n_estimators), random_state=0)
        _transition_model.fit(self.X, self.y)
        _transition_model_select = SelectFromModel(
            _transition_model, prefit=True, threshold=top_args
        )
        _transition_model_x = _transition_model_select.transform(self.X)

        # In this case there are not features to select, therefore, return 0 as accuracy.
        if _transition_model_x.shape[1] == 0:
            return 0

        # score = make_scorer(self.loss_function_2, greater_is_better=True)
        # af = AffinityPropagation().fit(self.X)
        # labels = af.labels_
        # val = metrics.metrics.adjusted_mutual_info_score(self.y, )
        # val = cross_validate(
        # AffinityPropagation(),
        # X=_transition_model_x,
        # y=self.y,
        # scoring=score,
        # cv=2
        # )

        score = make_scorer(self.loss_function_2, greater_is_better=True)
        val = cross_validate(
            ETC(n_estimators=int(n_estimators), random_state=0),
            X=_transition_model_x,
            y=self.y,
            scoring=score,
            cv=2
        )

        self.cumulative_objective_function.append({
            "score": val['test_score'].mean(),
            "n_estimators": n_estimators,
            "top_args": _transition_model_x.shape[1],
            "importance_cutoff": top_args
        })

        return val['test_score'].mean()

    def optimize(self):
        self.gp_params = {"alpha": 1e-5}
        self.etc_0 = BayesianOptimization(
            self.etc_ccv,
            {
                'n_estimators': (1000, 1000),
                'top_args': (self.min_importance, self.max_importance)
            }
        )

        self.etc_0.maximize(n_iter=self.epochs, **self.gp_params)

        print('selecting best performance parameters ...')
        selected_parameters = sorted(
            self.etc_0.res, key=lambda i: i['target'])[-1]

        self.forest = ETC(
            n_estimators=int(
                selected_parameters['params']['n_estimators']
            ),
            random_state=0
        )

        self.forest.fit(self.X, self.y)
        self._selected_features_model = SelectFromModel(
            self.forest, prefit=True, threshold=selected_parameters['params']['top_args']
        )

        self.parameters = pd.DataFrame({
            "score": [i['score'] for i in self.cumulative_objective_function],
            "n_estimators": [i['n_estimators'] for i in self.cumulative_objective_function],
            "top_args": [i['top_args'] for i in self.cumulative_objective_function],
            "importance_cutoff": [i['importance_cutoff'] for i in self.cumulative_objective_function]
        })
        # print(json.dumps(self.cumulative_objective_function, indent=4))

        self.x_t_selected = self._selected_features_model.transform(self.X)

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
        print('selecting discriminatory ARGs ...')
        print('saving into '+self.output_file+'.xlsx file ...')
        # writing in new excel file
        writer = pd.ExcelWriter(self.output_file+'.xlsx')
        self.X.to_excel(writer, 'cleaned_raw_input')
        self.x_selected.to_excel(writer, 'discriminatory_args')
        self.importances.to_excel(writer, 'ARGs_importance')
        self.parameters.to_excel(writer, 'ARGs_importance_cutoffs')
        writer.save()


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('--input-file', required=True, help='input excel file')
@click.option('--output-file', required=True, help='output file where to store the results')
# @click.option('--disc', default=50, help='top N discriminative ARGs (default 50)')
@click.option('--min-reads', default=1, help='minimum number of reads on each ARG (default 1)')
@click.option('--epochs', default=50, help='number of iterations the optimization algorithm run (default 50)')
@click.option('--max-importance', default=0.01, help='maximum importance for search space (default 0.01)')
@click.option('--min-importance', default=1e-5, help='minimum importance for search space (default 1e-5)')
# @click.option('--optimize', default=False, help='minimum number of reads on each ARG (default 1)')
def process(input_file, output_file, min_reads, epochs, max_importance, min_importance):
    """
    This program subtract the top N (50 default) discriminatory antibiotic resistance genes from a set of metagenomics samples.
    Hyperparameters of the supervised machine learning algorithm (extra tree classifier) are automatically tuned using the bayesian optimization.
    """

    if not input_file:
        print('\nUsage: extrarg --help\n')
        exit()

    extra_arg = ExtraARG(input_file, output_file,
                         min_reads, epochs, max_importance, min_importance)

    print('loading input datasets ...')
    extra_arg.load_data()

    print('performing optimization ...')
    extra_arg.optimize()

    print('building model with optimized parameters ...')
    extra_arg.discriminate()

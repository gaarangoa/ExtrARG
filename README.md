## Instalation (macOS/Linux)
To instal ExtraARGs algorighm just follow this instructions:

        git clone https://github.com/gaarangoa/ExtrARG.git
        cd ./ExtrARG/src/
        pip install .

Check the instalation by typing:

        extrarg -h

Run a test:

        extrarg --input-file ../test/PIRE_INFLUENT.xlsx --output-file ../test/temp

## Instalation (Windows)
### Requirements

        python 2.7 or python 3.x
        git
        pip

## Instalation

        git clone https://github.com/gaarangoa/ExtrARG.git
        cd ./ExtrARG/src/
        python -m pip install .

Check the instalation by typing:

        extrarg -h

Run a test

        python -m extrarg --input-file ../test/PIRE_INFLUENT.xlsx --output-file ../test/temp

## Usage

        extrarg --help # (macOS/linux)
        python -m extrarg --help # (Windows)
        Usage: extrarg [OPTIONS]

        This program subtract the top N (50 default) discriminatory antibiotic
        resistance genes from a set of metagenomics samples. Hyperparameters of
        the supervised machine learning algorithm (extra tree classifier) are
        automatically tuned using the bayesian optimization.

        Options:
        --input-file TEXT       input excel file
        --output-file TEXT      output file where to store the results
        --min-reads INTEGER     minimum number of reads on each ARG (default 1)
        --epochs INTEGER        number of iterations the optimization algorithm run
                                (default 10)
        --max-importance FLOAT  maximum importance for search space (default 0.01)
        --min-importance FLOAT  minimum importance for search space (default 1e-5)
        -h, --help              Show this message and exit.
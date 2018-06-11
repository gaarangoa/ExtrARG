from setuptools import setup

setup(
    name='extrarg',
    version='0.1',
    py_modules=['extrarg'],
    install_requires=[
        'Click',
        'numpy',
        'pandas',
        'sklearn',
        'bayesian-optimization'
    ],
    entry_points='''
        [console_scripts]
        extrarg=extrarg:process
    ''',
)
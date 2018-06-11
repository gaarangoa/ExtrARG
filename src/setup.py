from setuptools import setup

setup(
    name='extrarg',
    version='0.1',
    py_modules=['extrarg'],
    install_requires=[
        'Click=6.7',
        'numpy=1.14.3',
        'pandas=0.22.0',
        'xlrd',
        'sklearn=0.19.1',
        'bayesian-optimization'
    ],
    entry_points='''
        [console_scripts]
        extrarg=extrarg:process
    ''',
)
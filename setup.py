from setuptools import setup
setup(
    name='curate_adata',
    version='0.0.99',
    url='https://github.com/emdann/sfaira-data-curator',
    author='Emma Dann',
    license='MIT',
    entry_points={
        'console_scripts': [
            'curate_adata=curate_adata:main'
        ]
    },
    install_requires=[
        "scanpy >= 1.9.1",
        "sfaira >= 0.3.12"
    ]
)

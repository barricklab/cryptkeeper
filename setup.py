from setuptools import setup, find_packages

setup(
    name="cryptkeeper",
    version="0.1",
    packages=find_packages(),

    # Project uses reStructuredText, so ensure that the docutils get
    # installed or upgraded on the target machine
    install_requires=['ostir', 'rhotermpredict', 'promotercalculator', 'plotly'],

    # metadata to display on PyPI
    author="Barrick Lab",
    author_email="jbarrick@cm.utexas.edu",
    description="Package for predicting functional components of DNA sequences",
    keywords="DNA, RNA, Cryptkeeper, Barrick, RBS, Biology, Bio",
    url="https://github.com/barricklab/cryptkeeper",   # project home page, if any
    project_urls={
        "Bug Tracker": "https://github.com/barricklab/cryptkeeper/issues",
        "Source Code": "https://github.com/barricklab/cryptkeeper",
    },
    classifiers=[
        "License :: GNU General Public License v2.0"
    ],
    entry_points={
        'console_scripts' : [
          'cryptkeeper = cryptkeeper.cryptkeeper:main',
        ],
    }

    # could also include long_description, download_url, etc.
)

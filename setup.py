from setuptools import setup, find_packages

setup(
    name="cryptkeeper",
    version="1.0.0",
    packages=find_packages(include=["cryptkeeper"]),
    # Project uses reStructuredText, so ensure that the docutils get
    # installed or upgraded on the target machine
    install_requires=[
        "ostir",
        "rhotermpredict",
        "promotercalculator",
        "bokeh",
    ],
    # metadata to display on PyPI
    author="Barrick Lab",
    author_email="jbarrick@cm.utexas.edu",
    description="Package for predicting functional components of DNA sequences",
    keywords="DNA, RNA, CryptKeeper, Barrick, RBS, Biology, Bio",
    url="https://github.com/barricklab/cryptkeeper",  # project home page, if any
    project_urls={
        "Bug Tracker": "https://github.com/barricklab/cryptkeeper/issues",
        "Source Code": "https://github.com/barricklab/cryptkeeper",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Operating System :: OS Independent",
    ],
    entry_points={
        "console_scripts": [
            "cryptkeeper = cryptkeeper.cryptkeeper:main",
        ],
    },
    include_package_data=True,
    # could also include long_description, download_url, etc.
)

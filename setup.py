#!/usr/bin/env python3

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = []

setup_requirements = ['pytest-runner']

test_requirements = ['pytest>=3']

setup(
    author='Jonathon Vandezande',
    author_email='jevandezande@gmail.com',
    python_requires='>=3.8',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
    ],
    description='Generates clusters for X-ray computations following the model of Maganas et al. (DOI: 10.1039/c3cp50709b).',
    install_requires=requirements,
    license='MIT license',
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='spectra',
    name='spectra',
    packages=find_packages(include=['cluster', 'cluster.*']),
    scripts=[],
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/jevandezande/cluster',
    version='0.1.0',
    zip_safe=False,
)

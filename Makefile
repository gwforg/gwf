.PHONY: \
	help \
	init \
	test \
	lint \
	coverage \
	docs \
	clean \
	package-pypi \
	publish-pypi \
	install-conda \
	package-conda \
	publish-conda

export PATH := $(HOME)/.conda/bin:$(PATH)

help:
	@echo "Usage:"
	@echo "    make help        show this message"
	@echo "    make init        create and setup development environment"
	@echo "    make test        run unit tests"
	@echo "    make lint        run lint checks"
	@echo "    make test        build docs"
	@echo "    make clean       remove temporary files"
	@echo "    make package     build source distribution and wheel"
	@echo "    make publish     publish distributions to pypi"

init:
	pip install -r requirements.txt
	pip install -e . --no-deps

test:
	coverage run --source gwf -m pytest

lint:
	flake8 src/gwf

coverage:
	coverage report

docs:
	$(MAKE) -C docs html

clean:
	rm -rf docs/_build .gwfconf.json build/ dist/ .gwf .pytest_cache

# PyPI

package:
	rm -rf build dist .egg requests.egg-info
	python setup.py sdist bdist_wheel

publish:
	twine upload dist/*

# Conda

install-conda:
	wget -q https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
	chmod +x miniconda.sh
	./miniconda.sh -b -p "${HOME}/miniconda"
	rm miniconda.sh

	conda config --set always_yes true
	conda config --set anaconda_upload no
	conda config --add channels gwforg
	conda update --yes --quiet -n base conda
	conda install --quiet --yes conda-build=3.0.* anaconda-client=1.6.*

package-conda:
	conda build --python "${TRAVIS_PYTHON_VERSION}" --output-folder conda-bld/ conda/
	conda convert --platform all conda-bld/*/*.tar.bz2 -o conda-bld/

publish-conda:
	anaconda -t "${ANACONDA_TOKEN}" upload --force --no-progress --user gwforg conda-bld/*/*.tar.bz2
	rm -rf conda-bld/

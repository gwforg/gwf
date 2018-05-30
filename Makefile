.PHONY: \
	init \
	test \
	lint \
	coverage \
	coverage-ci \
	docs \
	clean \
	package-pypi \
	publish-pypi \
	install-conda \
	package-conda \
	publish-conda \
	publish

init:
	pip install pipenv --upgrade
	pipenv install --dev --skip-lock

test:
	pipenv run coverage run --source gwf setup.py test

lint:
	pipenv run flake8 src/gwf

coverage:
	pipenv run coverage report

coverage-ci:
	pip install coveralls
	coveralls

docs:
	pipenv run $(MAKE) -C docs html

clean:
	rm -rf docs/_build .gwfconf.json

package-pypi:
	pip install twine
	python setup.py sdist bdist_wheel

publish-pypi: package-pypi
	twine upload dist/*
	rm -rf build dist .egg requests.egg-info

install-conda:
	wget -q https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
	chmod +x miniconda.sh
	./miniconda.sh -b -p "${HOME}/miniconda"
	rm miniconda.sh

	export PATH="${HOME}/miniconda/bin:${PATH}" && \
	conda config --set always_yes true && \
	conda config --set anaconda_upload no && \
	conda config --add channels gwforg && \
	conda update --yes --quiet -n base conda && \
	conda install --quiet --yes conda-build=3.0.* anaconda-client=1.6.*

package-conda: install-conda
	export PATH="${HOME}/miniconda/bin:${PATH}" && \
	conda build --python "${TRAVIS_PYTHON_VERSION}" --output-folder conda-bld/ conda/ && \
	conda convert --platform all conda-bld/*/*.tar.bz2 -o conda-bld/

publish-conda: package-conda
	export PATH="${HOME}/miniconda/bin:${PATH}" && \
	anaconda -t "${ANACONDA_TOKEN}" --verbose --show-traceback upload --force --no-progress --user gwforg conda-bld/*/*.tar.bz2
	rm -rf conda-bld/

publish: publish-pypi publish-conda
	echo success

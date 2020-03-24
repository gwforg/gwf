.PHONY: \
	help \
	init \
	test \
	docs \
	clean \
	update-deps \
	init \
	update \
	package \
	publish \
	publish-test \
	package-conda \
	publish-conda

update-deps:
	pip install --upgrade pip-tools pip setuptools
	pip-compile --upgrade --build-isolation --generate-hashes --output-file requirements.txt
	pip-compile --upgrade --build-isolation --generate-hashes --output-file requirements-dev.txt requirements-dev.in

init:
	pip install --upgrade -r requirements.txt -r requirements-dev.txt
	pip install --no-deps --editable .
	rm -rf .nox

update: update-deps init

test:
	nox

clean:
	find . -prune -name ".egg-info" -type d -exec rm -rf {} ';'
	find . -prune -name ".gwf" -type d -exec rm -rf {} ';'
	find . -prune -name ".eggs" -type d -exec rm -rf {} ';'
	find . -prune -name "__pycache__" -type d -exec rm -rf {} ';'
	rm -rf docs/_build .gwfconf.json build/ dist/ .gwf .pytest_cache .egg conda-bld .coverage .nox

# PyPI

package: clean
	pip install --upgrade pep517
	rm -rf build/ dist/
	python -m pep517.build .

publish-test:
	twine upload --repository-url https://test.pypi.org/legacy/ dist/*

publish:
	twine upload --repository-url https://upload.pypi.org/legacy/ dist/*

# Conda

package-conda: clean
	conda build --output-folder conda-bld/ conda/

publish-conda:
	anaconda -t "${ANACONDA_TOKEN}" upload --force --no-progress --user gwforg conda-bld/*/*.tar.bz2

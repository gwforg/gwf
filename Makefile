.PHONY: \
	clean \
	build \
	publish \
	publish-test \
	package-conda \
	publish-conda

clean:
	find . -prune -name ".egg-info" -type d -exec rm -rf {} ';'
	find . -prune -name ".gwf" -type d -exec rm -rf {} ';'
	find . -prune -name ".eggs" -type d -exec rm -rf {} ';'
	find . -prune -name "__pycache__" -type d -exec rm -rf {} ';'
	rm -rf docs/_build .gwfconf.json build/ dist/ .gwf .pytest_cache .egg conda-bld .coverage .nox

build:
	flit build --no-setup-py --format wheel

# PyPI

publish-test:
	flit publish --repository testpypi dist/*

publish:
	flit publish --repository pypi dist/*

# Conda

package-conda: clean
	rm -rf conda-bld/
	conda build --output-folder conda-bld/ conda/

publish-conda:
	anaconda -t "${ANACONDA_TOKEN}" upload --force --no-progress --user gwforg conda-bld/*/*.tar.bz2

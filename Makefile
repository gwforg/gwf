.PHONY: init test lint coverage docs

init:
	pip install -r requirements.txt
	pip install -e .

test:
	coverage run --source gwf setup.py test

lint:
	flake8

coverage:
	coverage report

docs:
	$(MAKE) -C docs html

install-slurm:
	sudo apt-get update
	sudo apt-get install -y -qq slurm-llnl
	sudo cp ci/slurm.conf /etc/slurm-llnl/slurm.conf
	sudo mkdir -p /var/{run,spool,slurmd}
	sudo /usr/sbin/create-munge-key -f
	sudo /etc/init.d/munge restart
	sudo /etc/init.d/slurm-llnl restart

integration-test:
	cd examples/minimal-workflow/ && gwf

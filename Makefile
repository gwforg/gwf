.PHONY: init test lint coverage docs init-slurm

init:
	pip install -r requirements.txt
	pip install -e .

init-slurm:
	sudo apt-get update
	sudo apt-get install -y -qq slurm-llnl
	sudo mkdir -p /var/{run,spool,slurmd}
	sudo cp ci/slurm.conf /etc/slurm-llnl/
	sudo /usr/sbin/create-munge-key -f
	sudo /etc/init.d/munge restart
	sudo /etc/init.d/slurm-llnl restart

test:
	coverage run --source gwf setup.py test

lint:
	flake8

coverage:
	coverage report

docs:
	$(MAKE) -C docs html

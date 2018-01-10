init:
	pip install -r requirements.txt

test:
	nosetests -v --with-coverage --cover-tests tests
# Option Calculator Website

This repository contains sources of Django application that powers Option Calculator.

## What's in it?

It's a simple calculator that contains 7 models:

- European call/put option
- Implied volatility calculator
- American call/put option
- Geometric Asian option
- Arithmetic Asian option
- Geometric basket option
- Arithmetic basket option

## Vedio

<img src="https://cl.ly/452m130H1C1H/Screen%20Recording%202018-04-20%20at%2006.36%20%E4%B8%8B%E5%8D%88.gif" />



# Run the website

Please note that we use Python 2.7, so make sure that you use correct version when running commands below.

## Setting up a development environment

### Preparation for Windows

Download pip and install
go to https://pip.pypa.io/en/latest/installing/

Then install virtualenv and virtualenvwrapper

    pip install virtualenvwrapper-win
 
Then create a virtual environment for your project:

    mkvirtualenv myproject
    
If you start a new command prompt, you’ll need to activate the environment again using:

    workon myproject

### Steps for Windows/Mac/Ubuntu

First, clone the repository,
step into newly created `ass3_website` directory:

    cd ass3_website

Create a new virtual environment if needed. Then, install all the required dependencies:

    pip install -r requirements.txt


Run the migration

    python manage.py makemigrations
    python manage.py migrate

Run your local server:

    python manage.py runserver

If have any problem, you can try to edit the requirements.txt,or install/update the scipy/django/numpy manually.
Should you have any problem, you can contact Kexin Chen.

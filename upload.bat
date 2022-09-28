
@echo on
set path=D:\Anaconda3;D:\Anaconda3\Library\bin;D:\Anaconda3\Scripts;D:\Anaconda3\condabin;%path%
set path=C:\Users\Administrator\anaconda3;C:\Users\Administrator\anaconda3\Library\bin;%path%
set path=C:\Users\Administrator\anaconda3\Scripts;C:\Users\Administrator\anaconda3;\condabin;%path%


python setup.py sdist

twine check dist/*

twine upload dist/*

rd /s /Q dist

rd /s /Q ase_db_extension.egg-info

pause

pause

exit
require 5.005;
use inc::Module::Install;

name('Statistics-Regression');
author('Ivo Welch <ivo.welch@yale.edu>');
abstract('weighted linear regression package (line+plane fitting)');
license('GPL');
version_from('Regression.pm');

requires(qw(perl 5.005));
build_requires(qw(Test::More 0));
include_deps('Test::More', 5.005);
check_nmake();

&WriteAll();

FROM ubuntu:14.04

# Requires python 2 and 3
RUN apt-get update && apt-get -y install python3.4 python-pip python3-pip cpanminus
RUN apt-get update && apt-get -y install python-scipy python-lxml # cuz they don't work via pip
RUN pip install numpy openpyxl biopython flask
RUN pip3 install pysmb
RUN cpanm install Modern::Perl Getup::Long Data::Alias Filename::Basename Pod::Usage List::MoreUtils YAML::Syck Log::Log4perl File::NFSLock
CMD tail -f /dev/null

FROM ubuntu:14.04

# Requires python 2 and 3
RUN apt-get update && apt-get -y install python3.4 python-pip python3-pip cpanminus libncurses5-dev pigz default-jre #zlib
RUN apt-get update && apt-get -y install python-scipy python-lxml # cuz they don't work via pip
RUN pip install numpy openpyxl biopython
RUN pip3 install pysmb flask
RUN cpanm install Modern::Perl Getopt::Long Data::Alias File::Basename Pod::Usage List::MoreUtils YAML::Syck Log::Log4perl File::NFSLock

RUN mkdir -p /root/src
COPY . /root/src
WORKDIR /root/src

CMD tail -f /dev/null

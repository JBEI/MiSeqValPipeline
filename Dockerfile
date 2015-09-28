FROM ubuntu:14.04

VOLUME /root

# Requires python 2 and 3
RUN apt-get update
RUN apt-get -y install python3.4 python-pip python3-pip cpanminus \
 libncurses5-dev pigz default-jre gunicorn samtools bwa curl jq ruby1.9.3
RUN apt-get -y install python-scipy python-lxml # cuz they don't work via pip
RUN pip install numpy openpyxl biopython pysmb flask
RUN pip3 install pysmb flask
RUN cpanm install Modern::Perl Getopt::Long Data::Alias File::Basename Pod::Usage List::MoreUtils YAML YAML::Syck Log::Log4perl File::NFSLock

RUN mkdir -p /src
WORKDIR /src
VOLUME /src

EXPOSE 80
ENV SMB_USERNAME=$SMB_USERNAME
ENV SMB_PASSWORD=$SMB_PASSWORD

CMD cd /src/website && gunicorn -b :80 seqval:app

FROM ubuntu:16.04
MAINTAINER sminot@fredhutch.org

# Install prerequisites
RUN apt update && \
	apt-get install -y build-essential wget unzip python2.7 python-dev python-pip bats zlib1g-dev bzip2

# Use /share as the working directory
RUN mkdir /share
WORKDIR /share

# Add files
RUN mkdir /usr/vsearch
ADD . /usr/vsearch

# Get the binary from the latest release
RUN cd /usr/vsearch && \
	wget --no-check-certificate  https://github.com/torognes/vsearch/releases/download/v2.4.4/vsearch-2.4.4-linux-x86_64.tar.gz && \
	tar xzvf vsearch-2.4.4-linux-x86_64.tar.gz && \
	cd vsearch-2.4.4-linux-x86_64 && \
	ln -s /usr/vsearch/vsearch-2.4.4-linux-x86_64/bin/vsearch /usr/local/bin

# Run tests and then remove the folder
RUN bats /usr/vsearch/tests/ && rm -r /usr/vsearch/tests/

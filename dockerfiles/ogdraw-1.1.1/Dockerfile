FROM ubuntu:22.04
RUN apt update && apt install -y wget git \
    build-essential \
    perl \
    libbio-perl-perl \
    libimage-magick-perl \
    libpostscript-simple-perl \
    libxml-generator-perl \
    libxml-simple-perl \
    cpanminus \
    vim \
    unzip

RUN cpanm Bio::Perl
RUN cpanm Bio::Restriction::Analysis
WORKDIR /opt
RUN wget https://chlorobox.mpimp-golm.mpg.de/GeneMap-1.1.1.tar.gz && \
    wget https://chlorobox.mpimp-golm.mpg.de/OGDraw_Config_10.zip && \
    wget https://chlorobox.mpimp-golm.mpg.de/ogd_xml_plastid.xml && \
    wget https://chlorobox.mpimp-golm.mpg.de/ogd_xml_mitochondrion.xml
RUN tar -xvzf GeneMap-1.1.1.tar.gz && \
    unzip OGDraw_Config_10.zip
RUN rm GeneMap-1.1.1.tar.gz OGDraw_Config_10.zip
WORKDIR /opt/GeneMap-1.1.1
RUN perl Makefile.PL && \
    make && \
    make test && \
    make install
RUN sed -i 's:<policy domain="coder" rights="none" pattern="PS" />:<policy domain="coder" rights="read|write" pattern="PS" />:g' /etc/ImageMagick-6/policy.xml
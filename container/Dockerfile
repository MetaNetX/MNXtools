FROM ubuntu:20.04 AS builder

ARG MNXTOOLS_VERSION
ARG BUILD_DATE
ARG MNXREF_VERSION

# METADATA
# Use Open Containers Initiative (OCI)
# See https://rehansaeed.com/docker-labels-depth/
# See https://github.com/opencontainers/image-spec/blob/master/annotations.md
# Exist also another structured label initiative: http://label-schema.org/ superceded by OCI now
LABEL org.opencontainers.image.title="MNXtools: MetaNetX/MNXref utilities"
LABEL org.opencontainers.image.version=$MNXTOOLS_VERSION
LABEL org.opencontainers.image.vendor="MetaNetX.org"
LABEL org.opencontainers.image.authors="sebastien.moretti@sib.swiss"
LABEL org.opencontainers.image.url="https://www.metanetx.org/"
LABEL org.opencontainers.image.source="https://github.com/MetaNetX/MNXtools"
LABEL org.opencontainers.image.documentation="https://www.metanetx.org/"
LABEL org.opencontainers.image.licenses="GPL-3+"
LABEL org.opencontainers.image.description="MNXtools is a set of tools to use the MetaNetX/MNXref namespace \
as well as to manipulate models in MetaNetX.org"
LABEL org.opencontainers.image.created=$BUILD_DATE
LABEL mnx.mnxref.version=$MNXREF_VERSION
LABEL mnx.sbml.libsbml.version="5.19.0"
LABEL mnx.lang.perl.version="5.30.0"
LABEL mnx.lang.r.version="4.1.2"


# INSTALLATION
ENV USER=metanetx
# In a single layer to save time AND space
RUN echo '# OS update' \
	&& export DEBIAN_FRONTEND=noninteractive \
	&& apt-get update  -y \
	&& apt-get install -y --no-install-recommends perl libperl-dev libipc-run-perl libsort-naturally-perl libmodule-runtime-perl liblog-log4perl-perl libdigest-crc-perl libtest-simple-perl libmath-bigint-perl wget cmake make libexpat1-dev swig libyaml-dev libbz2-dev zlib1g-dev liblzma-dev check gcc g++ gfortran python3 libexpat1 libbz2-1.0 zlib1g liblzma5 parallel patchelf gnupg apt-utils ca-certificates libterm-readline-gnu-perl readline-common software-properties-common dirmngr git \
	&&        echo '# Install R (in a non-interactive way)' \
	&& echo '# See https://cran.r-project.org/bin/linux/ubuntu/' \
	&& echo '# Signature of Ubuntu archives on CRAN' \
	&& wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc \
	&& echo '# Add the R 4.1 repo from CRAN'  \
	&& add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/" \
	&& apt-get update -y \
	&& apt-get install -y  build-essential libssl-dev libncurses5-dev libcurl4-openssl-dev liblzma-dev libbz2-dev libboost-all-dev sqlite3 libsqlite3-0 libsqlite3-dev libgsl0-dev zlib1g-dev libxml2-dev libgmp-dev libgmp10 libmpfr-dev \
	&& apt-get install -y  r-base r-base-dev \
	&&        echo '# Install MNXtools' \
	&& cd /usr/local \
	&& git clone https://github.com/MetaNetX/MNXtools.git \
	&& cd MNXtools/bin/ \
	&& ln -s $PWD/* /usr/local/bin/ \
	&& cd .. \
	&& rm -Rf .dvc* .git* \
	&& cd cache \
	&& wget "https://www.metanetx.org/ftp/$MNXREF_VERSION/ChemSpace.bindump.gz" \
	&& gunzip ChemSpace.bindump.gz \
	&&        echo '# Install extra Perl lib' \
	&& cpan YAML::XS \
	&&        echo '# Install Perl libSBML' \
	&& cd /usr/local/ \
	&& wget 'https://downloads.sourceforge.net/project/sbml/libsbml/5.19.0/stable/libSBML-5.19.0-core-plus-packages-src.tar.gz' \
	&& tar xfz libSBML-*-core-plus-packages-src.tar.gz \
	&& ls -1 libSBML-*-core-plus-packages-src.tar.gz >>MANIFEST \
	&& rm -f libSBML-*-core-plus-packages-src.tar.gz \
	&& cd libSBML-*-Source/ \
	&& mkdir build \
	&& cd build/ \
	&& cmake -DWITH_JAVA=OFF -DWITH_CHECK=ON -DWITH_SWIG=ON -DWITH_PERL=ON -DWITH_PYTHON=OFF -DWITH_EXPAT=ON -DWITH_LIBXML=OFF -DENABLE_LAYOUT=ON -DENABLE_QUAL=ON -DENABLE_COMP=ON -DENABLE_FBC=ON -DENABLE_RENDER=ON -DENABLE_GROUPS=ON -DENABLE_MULTI=ON -DENABLE_DISTRIB=ON .. \
	&& make -j2 \
	&& make install \
	&& echo '## Fix for weird path that makes Perl LibSBML.pm unables to read LibSBML.so' \
	&& cd /usr/local/lib/perl5/site_perl/5.30.0/x86_64-linux-gnu-thread-multi/auto/ \
	&& ln -s libSBML LibSBML \
	&& echo '## Put in default Perl lib path' \
	&& ln -s $PWD/libSBML $PWD/LibSBML /usr/lib/x86_64-linux-gnu/perl/5.30.0/auto/ \
	&& cd .. \
	&& ln -s $PWD/LibSBML* /usr/lib/x86_64-linux-gnu/perl/5.30.0/ \
	&& cd /usr/local/libSBML-*-Source/build \
	&& make check \
	&& cd /usr/local \
	&& rm -Rf libSBML-*-Source/ \
	&& echo '# CLEANING' \
	&& echo '## So try to remove manually useless stuff' \
	&& apt remove -y *-dev wget cmake make swig autoconf check gcc g++ gfortran python3* patchelf \
	&& apt autoremove -y \
	&& apt-get clean \
	&& rm -rf /var/lib/apt/lists/* \
	&& groupadd -g 1001 ${USER} && useradd -r -u 1001 -g ${USER} ${USER}


# SECURITY
## Control root access
USER ${USER}
ENV HOME=/home/${USER}
WORKDIR ${HOME}
##Rest of Dockerfile with this user
ENV LC_ALL="C"
ENV MNXTOOLS_VERSION=$MNXTOOLS_VERSION
ENV MNXREF_VERSION=$MNXREF_VERSION


# TEST
RUN /usr/local/bin/SBML_to_TSV.pl || true

CMD ["/usr/local/bin/help"]


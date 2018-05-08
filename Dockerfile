FROM batk/bat-testenv:centos7-root6

COPY . /opt/src/bat

RUN cd /opt/src/bat \
    && tools/install-bat.sh /opt/bat

ENV \
    PATH="/opt/bat/bin:$PATH" \
    LD_LIBRARY_PATH="/opt/bat/lib:$LD_LIBRARY_PATH" \
    CPATH="/opt/bat/include:$CPATH" \
    PKG_CONFIG_PATH="/opt/bat/lib/pkgconfig:$PKG_CONFIG_PATH"

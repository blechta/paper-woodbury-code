FROM python:3.8.12-buster

# Install vtk build featuring XDMF3
WORKDIR /tmp
RUN wget -nv -O ParaView-5.9.1-osmesa-MPI-Linux-Python3.8-64bit.tar.gz "https://www.paraview.org/paraview-downloads/download.php?submit=Download&version=v5.9&type=binary&os=Linux&downloadFile=ParaView-5.9.1-osmesa-MPI-Linux-Python3.8-64bit.tar.gz" \
 && mkdir -p /usr/local/lib/paraview-5.9 \
 && tar -C /usr/local/lib/paraview-5.9 --strip-components=2 -xzf ParaView-5.9.1-osmesa-MPI-Linux-Python3.8-64bit.tar.gz ParaView-5.9.1-osmesa-MPI-Linux-Python3.8-64bit/lib/ \
 && find /tmp -mindepth 1 -delete
ENV PYTHONPATH="/usr/local/lib/paraview-5.9/python3.8/site-packages:${PYTHONPATH}"
ENV LD_LIBRARY_PATH="/usr/local/lib/paraview-5.9:${LD_LIBRARY_PATH}"

# Install pyvista
WORKDIR /tmp
RUN pip3 install -v --no-cache-dir --no-deps aiohttp==3.8.0 aiosignal==1.2.0 appdirs==1.4.4 async-timeout==4.0.0 attrs==21.2.0 charset-normalizer==2.0.7 cycler==0.11.0 frozenlist==1.2.0 idna==3.3 imageio==2.10.3 kiwisolver==1.3.2 matplotlib==3.4.3 meshio==4.4.6 multidict==5.2.0 pillow==8.4.0 pyparsing==3.0.4 python-dateutil==2.8.2 pyvista==0.32.1 scooby==0.5.7 six==1.16.0 typing-extensions==3.10.0.2 wslink==1.1.0 yarl==1.7.2 \
 && find /tmp -mindepth 1 -delete

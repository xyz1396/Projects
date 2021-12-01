TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

INCLUDEPATH += \
    src/MSToolkit/include/

HEADERS += \
    src/MSToolkit/include/MSNumpress.hpp \
    src/MSToolkit/include/MSObject.h \
    src/MSToolkit/include/MSReader.h \
    src/MSToolkit/include/MSToolkitTypes.h \
    src/MSToolkit/include/RAWReader.h \
    src/MSToolkit/include/Spectrum.h \
    src/MSToolkit/include/amigaconfig.h \
    src/MSToolkit/include/ascii.h \
    src/MSToolkit/include/asciitab.h \
    src/MSToolkit/include/crc32.h \
    src/MSToolkit/include/deflate.h \
    src/MSToolkit/include/expat.h \
    src/MSToolkit/include/expat_config.h \
    src/MSToolkit/include/expat_external.h \
    src/MSToolkit/include/gzguts.h \
    src/MSToolkit/include/iasciitab.h \
    src/MSToolkit/include/inffast.h \
    src/MSToolkit/include/inffixed.h \
    src/MSToolkit/include/inflate.h \
    src/MSToolkit/include/inftrees.h \
    src/MSToolkit/include/internal.h \
    src/MSToolkit/include/latin1tab.h \
    src/MSToolkit/include/macconfig.h \
    src/MSToolkit/include/mzMLWriter.h \
    src/MSToolkit/include/mzParser.h \
    src/MSToolkit/include/nametab.h \
    src/MSToolkit/include/pepXMLWriter.h \
    src/MSToolkit/include/sqlite3.h \
    src/MSToolkit/include/trees.h \
    src/MSToolkit/include/utf8tab.h \
    src/MSToolkit/include/winconfig.h \
    src/MSToolkit/include/xmlrole.h \
    src/MSToolkit/include/xmltok.h \
    src/MSToolkit/include/xmltok_impl.h \
    src/MSToolkit/include/zconf.h \
    src/MSToolkit/include/zlib.h \
    src/MSToolkit/include/zutil.h \
    src/Scores/CometSearchMod.h \
    src/Scores/MVH.h \
    src/SiprosReader.h \
    src/directoryStructure.h \
    src/isotopologue.h \
    src/ms2scan.h \
    src/ms2scanvector.h \
    src/peptide.h \
    src/proNovoConfig.h \
    src/proteindatabase.h \
    src/ptm.h \
    src/tokenvector.h

SOURCES += \
    src/MSToolkit/src/MSToolkit/MSObject.cpp \
    src/MSToolkit/src/MSToolkit/MSReader.cpp \
    src/MSToolkit/src/MSToolkit/RAWReader.cpp \
    src/MSToolkit/src/MSToolkit/Spectrum.cpp \
    src/MSToolkit/src/MSToolkit/mzMLWriter.cpp \
    src/MSToolkit/src/MSToolkit/pepXMLWriter.cpp \
    src/MSToolkit/src/expat-2.2.0/xmlparse.c \
    src/MSToolkit/src/expat-2.2.0/xmlrole.c \
    src/MSToolkit/src/expat-2.2.0/xmltok.c \
    src/MSToolkit/src/expat-2.2.0/xmltok_impl.c \
    src/MSToolkit/src/expat-2.2.0/xmltok_ns.c \
    src/MSToolkit/src/mzParser/BasicChromatogram.cpp \
    src/MSToolkit/src/mzParser/BasicSpectrum.cpp \
    src/MSToolkit/src/mzParser/Czran.cpp \
    src/MSToolkit/src/mzParser/MSNumpress.cpp \
    src/MSToolkit/src/mzParser/PWIZface.cpp \
    src/MSToolkit/src/mzParser/RAMPface.cpp \
    src/MSToolkit/src/mzParser/mz5handler.cpp \
    src/MSToolkit/src/mzParser/mzMLReader.cpp \
    src/MSToolkit/src/mzParser/mzParser.cpp \
    src/MSToolkit/src/mzParser/mzpMz5Config.cpp \
    src/MSToolkit/src/mzParser/mzpMz5Structs.cpp \
    src/MSToolkit/src/mzParser/mzp_base64.cpp \
    src/MSToolkit/src/mzParser/saxhandler.cpp \
    src/MSToolkit/src/mzParser/saxmzmlhandler.cpp \
    src/MSToolkit/src/mzParser/saxmzxmlhandler.cpp \
    src/MSToolkit/src/sqlite-3.7.7.1/sqlite3.c \
    src/MSToolkit/src/sqlite-3.7.7.1/sqlite3.cpp \
    src/MSToolkit/src/zlib-1.2.8/adler32.c \
    src/MSToolkit/src/zlib-1.2.8/compress.c \
    src/MSToolkit/src/zlib-1.2.8/crc32.c \
    src/MSToolkit/src/zlib-1.2.8/deflate.c \
    src/MSToolkit/src/zlib-1.2.8/gzclose.c \
    src/MSToolkit/src/zlib-1.2.8/gzlib.c \
    src/MSToolkit/src/zlib-1.2.8/gzread.c \
    src/MSToolkit/src/zlib-1.2.8/gzwrite.c \
    src/MSToolkit/src/zlib-1.2.8/infback.c \
    src/MSToolkit/src/zlib-1.2.8/inffast.c \
    src/MSToolkit/src/zlib-1.2.8/inflate.c \
    src/MSToolkit/src/zlib-1.2.8/inftrees.c \
    src/MSToolkit/src/zlib-1.2.8/trees.c \
    src/MSToolkit/src/zlib-1.2.8/uncompr.c \
    src/MSToolkit/src/zlib-1.2.8/zutil.c \
    src/Scores/CometSearchMod.cpp \
    src/Scores/MVH.cpp \
    src/SiprosReader.cpp \
    src/directoryStructure.cpp \
    src/isotopologue.cpp \
    src/main.cpp \
    src/mpimain.cpp \
    src/ms2scan.cpp \
    src/ms2scanvector.cpp \
    src/peptide.cpp \
    src/proNovoConfig.cpp \
    src/proteindatabase.cpp \
    src/ptm.cpp \
    src/tokenvector.cpp

DISTFILES += \
    src/MSToolkit/src/expat-2.2.0/COPYING \
    src/MSToolkit/src/zlib-1.2.8/ChangeLog \
    src/MSToolkit/src/zlib-1.2.8/README


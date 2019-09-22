
/** 
 * \class CDR3SeqData CDR3SeqData.h
 * \brief Class to store CDR3 information of a sequence.
 * \author C. Olivares
 * 
 */

#ifndef CDR3SEQDATA_H
#define CDR3SEQDATA_H

class CDR3SeqData {
public:
    CDR3SeqData();
    CDR3SeqData(const CDR3SeqData& orig);
    virtual ~CDR3SeqData();
    int seq_index;
    int v_anchor;
    int j_anchor;
    std::string CDR3nt;
    std::string CDR3aa;
    std::string strData();

private:

};

#endif /* CDR3SEQDATA_H */


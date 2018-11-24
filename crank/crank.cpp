#include "stdafx.h"
#include "crank.h"

int main(int argc, char* argv[]) {
  Env = TEnv(argc, argv, TNotify::StdNotify);
  Env.PrepArgs(TStr::Fmt("crank. build: %s, %s. Time: %s", __TIME__, __DATE__, TExeTm::GetCurTm()));
  TExeTm ExeTm;
  Try
  const TStr InFNa = Env.GetIfArgPrefixStr("-c:", "karate_communities.txt", "Community affiliation data");
  const TStr InFNm = Env.GetIfArgPrefixStr("-i:", "karate.txt", "Input edgelist file name");
  const TStr OutFNm = Env.GetIfArgPrefixStr("-o:", "prioritization.txt", "Output file name");
  const bool FitProba = Env.GetIfArgPrefixBool("-p:", true, "Fit an auxiliary network model (works with non-statistical community detection methods)");
  const TStr InFNn = Env.GetIfArgPrefixStr("-in:", "", "Input file name (probabilities of nodes belonging to communities)");
  const TStr InFNe = Env.GetIfArgPrefixStr("-ie:", "", "Input file name (probabilities of edges)");
  const TStr InFNc = Env.GetIfArgPrefixStr("-ic:", "", "Input file name (probabilities of edges given the communities)");
  const double Alpha = Env.GetIfArgPrefixFlt("-a:", 0.15, "Network perturbation intensity alpha");
  const double Prior = Env.GetIfArgPrefixFlt("-pr:", 0.05, "Relative size of temporary gold standard (p ~ prior probability)");
  int NBin = Env.GetIfArgPrefixInt("-b:", -1, "Number of bins (-1: detect automatically, B = |C|/50)");
  const int MaxIter = Env.GetIfArgPrefixInt("-mx:", 20, "Maximum number of iterations for rank aggregation");
  const double CorrStop = Env.GetIfArgPrefixFlt("-s:", 1., "Convergence criterion for rank aggregation");
  const TStr InFNgs = Env.GetIfArgPrefixStr("-gs:", "", "Input file name (side information)");
  TRnd Rnd(10);
  // Loading edgelist
  PUNGraph G;
  TStrHash<TInt> NodeNameH;
  G = TCRankUtil::LoadEdgeListStr<PUNGraph>(InFNm, NodeNameH);
  printf("Graph: %d Nodes %d Edges\n", G->GetNodes(), G->GetEdges());
  // Loading community affiliation data
  TVec<TIntV> CmtyVV;
  TVec<TStr> CmtyNames;
  THash<TStr, TInt> CmtyNameH;
  int Membs = 0, CId = 0;
  TSsParser Ss(InFNa, ssfWhiteSep);
  while (Ss.Next()) {
    if(Ss.GetFlds() > 0) {
      TIntV CmtyV;
      CmtyNames.Add(Ss.GetFld(0));
      CmtyNameH.AddDat(Ss.GetFld(0), CId);
      for(int i = 1; i < Ss.GetFlds(); i++) {
        char *Name = Ss[i];
        int NId = NodeNameH.GetKeyId(Name);
        CmtyV.Add(NId);
      }
      CmtyVV.Add(CmtyV);
      Membs += CmtyV.Len();
      CId++;
    }
  }
  printf("Community loading completed (%d communities, %d memberships)\n", CmtyVV.Len(), Membs);
  // Loading side (orthogonal, external) information
  // We assume that communities are given in the same order as in -c: file
  TFltV SideV;
  SideV.Gen(CmtyVV.Len());
  if (InFNgs != "") {
    TSsParser Sgs(InFNgs, ssfWhiteSep);
    while (Sgs.Next()) {
      if(Sgs.GetFlds() > 0) {
        int CId = CmtyNameH.GetDat(Sgs.GetFld(0));
        SideV[CId] = Sgs.GetFlt(1);
        printf("%s\t%d\t%f\n", Sgs.GetFld(0), CId, SideV[CId].Val);
      }
    }
  }
  printf("Side information loading completed\n");
  // Loading edge and community affiliation probabilities (if available)
  THash<TIntPr,TFlt> EdgeProbaH;
  THash<TIntPr,TFlt> CProbaH;
  THash<TIntTr,TFlt> CEdgeProbaH;
  if (! FitProba) {
    // edge probabilities
    TSsParser Sse(InFNe, ssfWhiteSep);
    while (Sse.Next()) {
      if(Sse.GetFlds() > 0) {
        char *SrcName = Sse[0];
        int SrcNId = NodeNameH.GetKeyId(SrcName);
        char *DstName = Sse[1];
        int DstNId = NodeNameH.GetKeyId(DstName);
        float Proba = Sse.GetFlt(2);
        EdgeProbaH.AddDat(TIntPr(SrcNId, DstNId), Proba);
      }
    }
    printf("Edge probability loading completed (%d probabilities)\n", EdgeProbaH.Len());
    // probabilities of edges given the communities
    TSsParser Ssc(InFNc, ssfWhiteSep);
    while (Ssc.Next()) {
      if(Ssc.GetFlds() > 0) {
        char *SrcName = Ssc[0];
        int SrcNId = NodeNameH.GetKeyId(SrcName);
        char *DstName = Ssc[1];
        int DstNId = NodeNameH.GetKeyId(DstName);
        if (!CmtyNameH.IsKey(Ssc.GetFld(2))) {
          printf("Community: %s has no membership data\n", Ssc.GetFld(2));
        }
        else {
          int CId = CmtyNameH.GetDat(Ssc.GetFld(2));
          float Proba = Ssc.GetFlt(3);
          CEdgeProbaH.AddDat(TIntTr(SrcNId, DstNId, CId), Proba);
        }
      }
    }
    printf("Community-specific edge probability loading completed (%d probabilities)\n", CEdgeProbaH.Len());
    // probabilities of nodes belonging to communities
    TSsParser Ssn(InFNn, ssfWhiteSep);
    while (Ssn.Next()) {
      if(Ssn.GetFlds() > 0) {
        char *SrcName = Ssn[0];
        int SrcNId = NodeNameH.GetKeyId(SrcName);
        if (!CmtyNameH.IsKey(Ssn.GetFld(1))) {
          printf("Community: %s has no membership data\n", Ssn.GetFld(1));
        }
        else {
          int CId = CmtyNameH.GetDat(Ssn.GetFld(1));
          float Proba = Ssn.GetFlt(2);
          CProbaH.AddDat(TIntPr(SrcNId, CId), Proba);
        }
      }
    }
    printf("Node-community probability loading completed (%d probabilities)\n", CProbaH.Len());
  }
  // Ranking communities
  TFltV RndV;
  TFltV PhiV;
  TFltV ModV;
  TFltV CRankV;
  TCRankUtil::GetRandomCmty(RndV, G, CmtyVV, Rnd);
  TCRankUtil::GetConductanceCmty(PhiV, G, CmtyVV);
  TCRankUtil::GetModularityCmty(ModV, G, CmtyVV);
  TCRank::Prioritize(CRankV, G, CmtyVV, Rnd, Alpha, Prior, NBin, MaxIter, CorrStop, FitProba, CProbaH, EdgeProbaH, CEdgeProbaH, SideV);
  // Saving results
  FILE* F = fopen(OutFNm.CStr(), "wt");
  fprintf(F, "Community\tCRank\tConductance\tModularity\tRandom\n");
  for (int i = 0; i < CmtyVV.Len(); i++) {
    fprintf(F, "%s\t%f\t%f\t%f\t%f\n", CmtyNames[i].CStr(), CRankV[i].Val, PhiV[i].Val, ModV[i].Val, RndV[i].Val);
  }
  fclose(F);
  Catch
  printf("run time: %s (%s)\n", ExeTm.GetTmStr(), TSecTm::GetCurTm().GetTmStr().CStr());
}

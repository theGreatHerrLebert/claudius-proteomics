import axios from 'axios';

// Use direct backend URL for static file serving, or /api for dev server
const API_BASE = import.meta.env.DEV ? '/api' : 'http://localhost:8000';

const api = axios.create({
  baseURL: API_BASE,
});

export interface PrecursorSummary {
  precursor_id: number;
  raw_file: string;
  mz: number;
  charge: number;
  rt_seconds: number;
  mobility: number;
  n_engines: number;
  confidence_weight: number | null;
  // FragPipe
  fragpipe_peptide: string | null;
  fragpipe_modified: string | null;
  fragpipe_probability: number | null;
  fragpipe_qvalue: number | null;
  // Sage
  sage_peptide: string | null;
  sage_modified: string | null;
  sage_qvalue: number | null;
  sage_match_tier: string | null;
  // DIA-NN
  diann_peptide: string | null;
  diann_modified: string | null;
  diann_qvalue: number | null;
  diann_ccs: number | null;
  diann_match_tier: string | null;
  diann_match_score: number | null;
  // Raw
  raw_intensity_meta: number | null;
  frame_id: number | null;
  isolation_mz: number | null;
  // Quality metrics
  ms1_rt_sigma: number | null;
  ms1_rt_r2: number | null;
  ms1_im_sigma: number | null;
  ms1_im_r2: number | null;
  isotope_cosim: number | null;
}

export interface RawFileInfo {
  name: string;
  count: number;
}

export interface PrecursorDetail extends PrecursorSummary {
  fragment_mz: number[];
  fragment_intensity: number[];
  fragment_mobility: number[];
  fragment_scan: number[];
  xic_rt: number[];
  xic_intensity: number[];
  mobilogram_im: number[];
  mobilogram_intensity: number[];
  isotope_mz: number[];
  isotope_intensity: number[];
  raw_rt: number[];
  raw_mz: number[];
  raw_mobility: number[];
  raw_intensity: number[];
}

export interface StoreInfo {
  path: string;
  num_precursors: number;
  num_row_groups: number;
  columns: string[];
}

export interface StoreStats {
  total_precursors: number;
  by_engines: Record<string, number>;
  by_charge: Record<string, number>;
  mz_range: [number, number];
  raw_files: Record<string, number>;
}

export async function loadStore(path: string) {
  const response = await api.post('/load', null, { params: { path } });
  return response.data;
}

export async function getStoreInfo(): Promise<StoreInfo> {
  const response = await api.get('/info');
  return response.data;
}

export async function getStats(): Promise<StoreStats> {
  const response = await api.get('/stats');
  return response.data;
}

export async function listPrecursors(params: {
  offset?: number;
  limit?: number;
  min_engines?: number;
  charge?: number;
  raw_file?: string;
  has_ms1?: boolean;
  sort_by?: string;
  sort_desc?: boolean;
}): Promise<PrecursorSummary[]> {
  const response = await api.get('/precursors', { params });
  return response.data;
}

export async function getRawFiles(): Promise<RawFileInfo[]> {
  const response = await api.get('/raw_files');
  return response.data;
}

export async function getPrecursor(id: number): Promise<PrecursorDetail> {
  const response = await api.get(`/precursor/${id}`);
  return response.data;
}

export default api;

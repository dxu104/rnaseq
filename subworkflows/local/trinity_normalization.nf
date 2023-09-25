//customized modules
//include { TRINITY as TRINITY_NORMALIZATION} from '../../modules/nf-core/Trinity/main'
include { TrinityNormalizeReads as TrinityNormalizeReads_SingleEnd } from '../../modules/local/TrinityNormalization/trinity_normalization.nf'
include { TrinityNormalizeReads as TrinityNormalizeReads_DoubleEnd } from '../../modules/local/TrinityNormalization/trinity_normalization.nf'


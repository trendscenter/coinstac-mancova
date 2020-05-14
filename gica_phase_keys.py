import coinstac_node_ops.local as ops_local
import coinstac_node_ops.remote as ops_remote
import coinstac_masking.local as mask_local
import coinstac_masking.remote as mask_remote
import coinstac_decentralized_row_means.local as drm_local
import coinstac_decentralized_row_means.remote as drm_remote
import coinstac_spatially_constrained_ica.local as scica_local
import coinstac_spatially_constrained_ica.remote as scica_remote
#import coinstac_backreconstruction.local as br_local
#import coinstac_backreconstruction.remote as br_remote
#import coinstac_decentralized_pca.local as dpca_local
#import coinstac_decentralized_pca.remote as dpca_remote
#import coinstac_gica.local as gica_local
#import coinstac_gica.remote as gica_remote
#import coinstac_ddfnc_preproc.local as dfncpp_local
#import coinstac_ddfnc_preproc.remote as dfncpp_remote
#import coinstac_dkmeans_ms.local as dkm_local
#import coinstac_dkmeans_ms.remote as dkm_remote

# Init

NOOP_LOCAL = [
    dict(
        do=[
            ops_local.local_noop
        ],
        recv=[],
        send='local_noop',
        args=[
            []
        ],
        kwargs=[
            {}
        ]
    )
]

NOOP_REMOTE = [
    dict(
        do=[
            ops_remote.remote_noop
        ],
        recv="local_noop",
        send='remote_noop',
        args=[
            []
        ],
        kwargs=[
            {}
        ]
    )
]

INIT_LOCAL = [
    dict(
        do=[
            ops_local.local_load_datasets,
            ops_local.local_output_to_input,
            mask_local.masking_local_1,
            ops_local.local_output_to_cache,
            ops_local.local_cache_to_input,
            drm_local.drm_local_1,
        ],
        recv=[],
        send='local_init',
        args=[
            [],
            [],
            [],
            [],
            []
        ],
        kwargs=[
            {},
            {},
            {},
            {},
            {}
        ],
    )
]

INIT_REMOTE = [
    dict(
        do=[
            drm_remote.drm_remote_1,
            ops_remote.remote_output_to_cache,
            ops_remote.remote_dump_cache_to_mat
        ],
        recv=INIT_LOCAL[0].get("send"),
        send='remote_init',
        args=[
            [],
            [],
            []
        ],
        kwargs=[
            {},
            {},
            {}
        ],

    )
]

# Spatially Constrained ICA
SPATIALLY_CONSTRAINED_ICA_LOCAL = [
    dict(
        do=[
            scica_local.scica_local_1,
            ops_local.local_output_to_cache
        ],
        recv=[],
        send='scica_local_1',
        args=[
            []
        ],
        kwargs=[
            {}
        ],
    )
]
SPATIALLY_CONSTRAINED_ICA_REMOTE = [
    dict(
        do=[
            ops_remote.remote_noop
        ],
        recv=SPATIALLY_CONSTRAINED_ICA_LOCAL[0].get('send'),
        send='scica_remote_noop',
        args=[
            []
        ],
        kwargs=[
            {}
        ],
    )
]

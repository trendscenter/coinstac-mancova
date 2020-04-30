import coinstac_node_ops.local as ops_local
import coinstac_node_ops.remote as ops_remote
import coinstac_mancova.local as mancova_local
import coinstac_mancova.remote as mancova_remote

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

MANCOVA_LOCAL = [
    dict(
        do=[
            # ops_local.local_load_datasets,
            # ops_local.local_output_to_input,
            mancova_local.local_run_mancova,
        ],
        recv=[],
        send='local_mancova_1',
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

MANCOVA_REMOTE = [
    dict(
        do=[
            ops_remote.remote_noop,
            # mancova_remote.mancova_aggregate
        ],
        recv=MANCOVA_LOCAL[0].get("send"),
        send='remote_aggregate',
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

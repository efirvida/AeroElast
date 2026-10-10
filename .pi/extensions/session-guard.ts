/**
 * session-guard — passive session-budget guidance for this repository.
 *
 * Why this exists (measured 2026-10-09 over 13 Pi sessions in this repo): 94 % of the
 * billed tokens were accumulated conversation history, tool payload was 0.15 %, and the
 * two marathon sessions (944 and 570 turns) were 64 % of the total. Cost grows with the
 * square of the session length, so the only cheap fix is closing at the right moment.
 *
 * This extension makes that moment visible without spending a single token: the hint lives
 * in the status bar plus at most three one-off notifications, never in the model context.
 *
 *   /handoff [feature]   ask the agent to leave the session handoff in odd/tasks/
 *
 * Project-local by design (`.pi/extensions/`), so no other repository is affected.
 * See odd/tasks/token-efficiency.md for the numbers behind the thresholds.
 */
import { existsSync, readdirSync } from "node:fs";

import type { ExtensionAPI, ExtensionContext } from "@earendil-works/pi-coding-agent";

const STATUS_KEY = "session-budget";

/** Ordered from the most severe; the first level whose threshold is reached wins. */
const LEVELS = [
    {
        at: 500_000,
        key: "marathon",
        label: "CARO: CERRÁ YA",
        type: "error" as const,
        notify:
            "Contexto >500k: sesión maratón. Cerrá ahora con /handoff — cada turno extra " +
            "re-manda todo el historial y es lo que se lleva la mayor parte del gasto.",
    },
    {
        at: 250_000,
        key: "close",
        label: "CERRÁ + /handoff",
        type: "warning" as const,
        notify:
            "Contexto >250k: cerrá en el próximo límite de work-unit. /handoff escribe el " +
            "relevo en odd/tasks/ y te da la receta para retomar en una sesión limpia.",
    },
    {
        at: 150_000,
        key: "watch",
        label: "cerrá pronto",
        type: "info" as const,
        notify:
            "Contexto >150k: empezá a pensar el cierre. Una sola cosa por sesión; el " +
            "relevo va en odd/tasks/<feature>.md.",
    },
];

/** The status-bar label for a context size, or "ok" below the first threshold. */
export function budgetLabel(tokens: number): string {
    for (const level of LEVELS) {
        if (tokens >= level.at) return level.label;
    }
    return "ok";
}

/** Compact status text: `ctx 132k · 41t · ok`. */
export function budgetStatus(tokens: number, turns: number): string {
    return `ctx ${Math.round(tokens / 1000)}k · ${turns}t · ${budgetLabel(tokens)}`;
}

export default function sessionGuard(pi: ExtensionAPI) {
    const announced = new Set<string>();
    let turns = 0;

    const render = (ctx: ExtensionContext, announce: boolean): void => {
        // Never let UI work break a turn: this extension is advisory only.
        try {
            const usage = ctx.getContextUsage();
            const tokens = usage?.tokens ?? null;
            if (tokens === null) return;

            if (ctx.hasUI) ctx.ui.setStatus(STATUS_KEY, budgetStatus(tokens, turns));

            if (!announce) return;
            const level = LEVELS.find((candidate) => tokens >= candidate.at);
            if (!level || announced.has(level.key)) return;
            announced.add(level.key);
            if (ctx.hasUI) ctx.ui.notify(level.notify, level.type);
        } catch {
            /* advisory extension: stay silent on failure */
        }
    };

    pi.on("session_start", (_event, ctx) => {
        announced.clear();
        turns = 0;
        render(ctx, false);
    });

    pi.on("turn_end", (_event, ctx) => {
        turns += 1;
        render(ctx, true);
    });

    pi.registerCommand("handoff", {
        description: "Cerrar la sesión dejando el relevo en odd/tasks/ (evita quemar tokens)",
        getArgumentCompletions: (prefix: string) => {
            // Cheap and useful: suggest the existing feature documents.
            try {
                if (!existsSync("odd/tasks")) return null;
                const names = readdirSync("odd/tasks").flatMap((name: string) =>
                    name.endsWith(".md") ? [name.replace(/\.md$/, "")] : [],
                );
                const matches = names.filter((name: string) => name.startsWith(prefix));
                return matches.length > 0 ? matches.map((value: string) => ({ value, label: value })) : null;
            } catch {
                return null;
            }
        },
        handler: async (args: string) => {
            const feature = args.trim();
            const target = feature
                ? `odd/tasks/${feature.replace(/\.md$/, "")}.md`
                : "odd/tasks/<feature>.md (usá el nombre de la feature activa, o creá el documento si no existe)";

            pi.sendUserMessage(
                [
                    "Cerrá esta sesión dejando el relevo, sin explorar nada nuevo:",
                    `1. Escribí o actualizá ${target} con: objetivo, estado actual, decisiones tomadas y por qué, evidencia verificada (comando + resultado), qué quedó a medias y el próximo paso concreto.`,
                    "2. Guardá en Engram sólo lo durable que aporte algo nuevo (topic_key + esencia), sin repetir el documento.",
                    "3. Respondé en tres líneas cómo retomar: qué leer, qué comando correr y cuál es el próximo paso.",
                    "No abras archivos nuevos, no corras tests y no repitas el estado de la conversación: usá lo que ya está en contexto.",
                ].join("\n"),
            );
        },
    });
}
